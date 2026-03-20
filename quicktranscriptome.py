#!/usr/bin/env python3
import argparse
import gzip
import math
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.request import urlretrieve

import pandas as pd


SPECIES_DEFAULTS = {
    "parapsilosis": {
        "fasta_urls": [
            "https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-60/fasta/candida_parapsilosis/dna/Candida_parapsilosis.GCA000182765v2.dna.toplevel.fa.gz",
            "https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-60/fasta/candida_parapsilosis/dna/Candida_parapsilosis.C_parapsilosis_CDC317.dna.toplevel.fa.gz",
            "http://www.candidagenome.org/download/sequence/C_parapsilosis_CDC317/C_parapsilosis_CDC317_current_chromosomes.fasta",
        ],
        "gff_urls": [
            "https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-60/gff3/candida_parapsilosis/Candida_parapsilosis.GCA000182765v2.60.gff3.gz",
            "https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-60/gff3/candida_parapsilosis/Candida_parapsilosis.C_parapsilosis_CDC317.60.gff3.gz",
            "http://www.candidagenome.org/download/gff/C_parapsilosis_CDC317/C_parapsilosis_CDC317_current_features.gff",
        ],
    }
}


def run_cmd(cmd):
    print("+", " ".join(cmd))
    subprocess.run(cmd, check=True)


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def should_skip(path: Path, resume: bool, label: str) -> bool:
    if resume and path.exists():
        print(f"Resume enabled: using existing {label}: {path}")
        return True
    return False


def download_file(url: str, out_path: Path):
    if out_path.exists():
        print(f"Using cached file: {out_path}")
        return
    print(f"Downloading {url} -> {out_path}")
    urlretrieve(url, out_path)


def normalize_urls(url_or_urls):
    if url_or_urls is None:
        return []
    if isinstance(url_or_urls, str):
        return [url_or_urls]
    return list(url_or_urls)


def download_first_available(urls, refs_dir: Path):
    errors = []
    for url in normalize_urls(urls):
        out_path = refs_dir / Path(url).name
        if out_path.exists():
            print(f"Using cached file: {out_path}")
            return out_path
        try:
            download_file(url, out_path)
            return out_path
        except (HTTPError, URLError) as exc:
            errors.append(f"{url} -> {exc}")
            print(f"Warning: could not download from {url} ({exc})")
    joined = "\n".join(errors) if errors else "No URLs provided."
    raise RuntimeError(f"Failed to download required reference file.\n{joined}")


def gunzip_if_needed(path: Path) -> Path:
    if path.suffix != ".gz":
        return path
    out = path.with_suffix("")
    if out.exists():
        return out
    with gzip.open(path, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out


def detect_annotation_format(annotation_path: Path) -> str:
    lower = annotation_path.name.lower()
    if lower.endswith(".gff") or lower.endswith(".gff3"):
        return "GFF"
    if lower.endswith(".gtf"):
        return "GTF"

    with open(annotation_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            attrs = parts[8]
            if "=" in attrs and "gene_id" not in attrs:
                return "GFF"
            return "GTF"
    return "GTF"


def featurecounts_annotation_args(annotation_path: Path):
    fmt = detect_annotation_format(annotation_path)
    if fmt == "GFF":
        # For many GFF3s, exon Parent points to transcript IDs; counting on gene features is safer.
        # Prefer gene_id when present (clean IDs), fallback to ID.
        gene_attr = "ID"
        with open(annotation_path, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                if parts[2] != "gene":
                    continue
                attrs = parts[8]
                if "gene_id=" in attrs:
                    gene_attr = "gene_id"
                break
        return ["-F", "GFF", "-t", "gene", "-g", gene_attr]
    return ["-F", "GTF", "-t", "exon", "-g", "gene_id"]


def canonical_sample_id(name: str) -> str:
    sample = str(name).strip()
    sample = Path(sample).name
    for suffix in [".sorted.bam", ".bam", ".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if sample.endswith(suffix):
            sample = sample[: -len(suffix)]
    for tag in ["_R1", "_R2"]:
        if tag in sample:
            sample = sample.split(tag)[0]
    return sample


def discover_samples(reads_dir: Path):
    reads = sorted(reads_dir.glob("*.fastq.gz")) + sorted(reads_dir.glob("*.fq.gz"))
    if not reads:
        raise FileNotFoundError(f"No FASTQ files found in {reads_dir}")

    paired = {}
    single = {}
    for fq in reads:
        name = fq.name
        if "_R1" in name:
            sample = name.split("_R1")[0]
            paired.setdefault(sample, {})["R1"] = fq
        elif "_R2" in name:
            sample = name.split("_R2")[0]
            paired.setdefault(sample, {})["R2"] = fq
        else:
            sample = fq.name.split(".")[0]
            single[sample] = fq

    samples = []
    for sample, mates in paired.items():
        if "R1" in mates and "R2" in mates:
            samples.append((sample, mates["R1"], mates["R2"]))
    for sample, fq in single.items():
        samples.append((sample, fq, None))

    if not samples:
        raise ValueError("FASTQ files found, but no valid samples detected.")
    return samples


def build_reference(species, refs_dir: Path, fasta_url=None, gff_url=None):
    ensure_dir(refs_dir)
    fasta_urls = normalize_urls(fasta_url)
    gff_urls = normalize_urls(gff_url)
    if species in SPECIES_DEFAULTS:
        default = SPECIES_DEFAULTS[species]
        if not fasta_urls:
            fasta_urls = default["fasta_urls"]
        if not gff_urls:
            gff_urls = default["gff_urls"]
    if not fasta_urls or not gff_urls:
        raise ValueError("You must provide --fasta-url and --gff-url for unknown species.")

    fasta_gz = download_first_available(fasta_urls, refs_dir)
    gff_gz = download_first_available(gff_urls, refs_dir)

    fasta = gunzip_if_needed(fasta_gz)
    gff = gunzip_if_needed(gff_gz)
    index = refs_dir / "reference.mmi"
    if not index.exists():
        run_cmd(["minimap2", "-d", str(index), str(fasta)])
    return fasta, gff, index


def align_and_count(args):
    out_dir = Path(args.out_dir)
    refs_dir = out_dir / "refs"
    bam_dir = out_dir / "bam"
    counts_dir = out_dir / "counts"
    ensure_dir(out_dir)
    ensure_dir(bam_dir)
    ensure_dir(counts_dir)

    _, gff, index = build_reference(
        species=args.species,
        refs_dir=refs_dir,
        fasta_url=args.fasta_url,
        gff_url=args.gff_url,
    )
    samples = discover_samples(Path(args.reads_dir))
    bam_files = []

    for sample, r1, r2 in samples:
        bam = bam_dir / f"{sample}.sorted.bam"
        bam_files.append(bam)
        bai = bam.with_suffix(bam.suffix + ".bai")
        if should_skip(bam, args.resume, "BAM"):
            if not bai.exists():
                run_cmd(["samtools", "index", str(bam)])
            continue

        if r2 is None:
            align_cmd = [
                "minimap2",
                "-ax",
                "splice",
                "-t",
                str(args.threads),
                str(index),
                str(r1),
            ]
        else:
            align_cmd = [
                "minimap2",
                "-ax",
                "splice",
                "-t",
                str(args.threads),
                str(index),
                str(r1),
                str(r2),
            ]

        samtools_sort = [
            "samtools",
            "sort",
            "-@",
            str(args.threads),
            "-o",
            str(bam),
            "-",
        ]
        print("+", " ".join(align_cmd), "|", " ".join(samtools_sort))
        p1 = subprocess.Popen(align_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(samtools_sort, stdin=p1.stdout)
        p1.stdout.close()
        rc2 = p2.wait()
        rc1 = p1.wait()
        if rc1 != 0 or rc2 != 0:
            raise subprocess.CalledProcessError(rc1 or rc2, align_cmd)
        run_cmd(["samtools", "index", str(bam)])

    count_file = counts_dir / "featureCounts.txt"
    if not should_skip(count_file, args.resume, "featureCounts output"):
        annotation_args = featurecounts_annotation_args(gff)
        run_cmd(
            [
                "featureCounts",
                "-T",
                str(args.threads),
                *annotation_args,
                "-a",
                str(gff),
                "-o",
                str(count_file),
                *[str(b) for b in bam_files],
            ]
        )

    # Convert featureCounts output to compact matrix.
    df = pd.read_csv(count_file, sep="\t", comment="#")
    fixed_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    sample_cols = [c for c in df.columns if c not in fixed_cols]
    clean = df[["Geneid", *sample_cols]].copy()
    clean.columns = ["gene_id", *[Path(c).stem.replace(".sorted", "") for c in sample_cols]]
    clean["gene_id"] = clean["gene_id"].astype(str).str.replace("^gene:", "", regex=True)
    clean_out = counts_dir / "counts_matrix.tsv"
    if not should_skip(clean_out, args.resume, "counts matrix"):
        clean.to_csv(clean_out, sep="\t", index=False)
        print(f"Counts matrix written: {clean_out}")

    counts_for_plots = clean.set_index("gene_id").T.astype(int)
    if args.metadata and args.condition_column:
        run_deseq(
            clean_out,
            Path(args.metadata),
            args.condition_column,
            counts_dir,
            args.resume,
            args.gsea_gmt,
            args.gsea_min_size,
            args.gsea_max_size,
        )
    else:
        write_pca(counts_for_plots, None, None, counts_dir, args.resume)
        write_heatmap(counts_for_plots, counts_dir, args.resume)


def _sample_ordered_counts_and_metadata(counts_matrix: Path, metadata_path: Path, condition_col: str):
    counts_df = pd.read_csv(counts_matrix, sep="\t")
    meta_df = pd.read_csv(metadata_path, sep=None, engine="python")
    meta_df.columns = [str(c).replace("\ufeff", "").strip() for c in meta_df.columns]
    meta_df = meta_df.loc[:, [c for c in meta_df.columns if not c.lower().startswith("unnamed:")]]

    if "sample" not in meta_df.columns:
        sample_candidates = [c for c in meta_df.columns if c.lower().strip() == "sample"]
        if sample_candidates:
            meta_df = meta_df.rename(columns={sample_candidates[0]: "sample"})
    if "sample" not in meta_df.columns:
        raise ValueError("Metadata must include a 'sample' column.")

    condition_col_clean = str(condition_col).strip()
    if condition_col_clean not in meta_df.columns:
        alt = [c for c in meta_df.columns if c.lower().strip() == condition_col_clean.lower()]
        if alt:
            condition_col_clean = alt[0]
        else:
            raise ValueError(f"Metadata missing condition column: {condition_col}")

    counts_df = counts_df.set_index("gene_id").T
    counts_df.index = pd.Index([canonical_sample_id(s) for s in counts_df.index], name="sample")
    counts_df = counts_df.astype(int)

    meta_df = meta_df.copy()
    meta_df["sample"] = meta_df["sample"].astype(str).str.strip()
    meta_df[condition_col_clean] = meta_df[condition_col_clean].astype("string").str.strip()
    meta_df["_sample_canonical"] = meta_df["sample"].map(canonical_sample_id)

    # Keep first metadata row per canonical sample to avoid ambiguous joins.
    meta_df = meta_df.drop_duplicates(subset="_sample_canonical", keep="first")
    meta_df = meta_df.set_index("_sample_canonical")

    common = counts_df.index.intersection(meta_df.index)
    if len(common) < 2:
        raise ValueError(
            "No overlapping samples between counts and metadata after normalization. "
            "Ensure metadata 'sample' values match FASTQ/BAM-derived sample names."
        )
    counts_df = counts_df.loc[common]
    meta_df = meta_df.loc[common]
    if meta_df[condition_col_clean].isna().all() or (meta_df[condition_col_clean] == "").all():
        fallback_cols = [c for c in ["treatment", "condition", "group"] if c in meta_df.columns and c != condition_col_clean]
        for fb in fallback_cols:
            if not (meta_df[fb].isna().all() or (meta_df[fb].astype("string").str.strip() == "").all()):
                print(
                    f"Warning: metadata column '{condition_col_clean}' is empty for matched samples; "
                    f"using '{fb}' instead."
                )
                condition_col_clean = fb
                break
    meta_df.index.name = "sample"
    return counts_df, meta_df, condition_col_clean


def print_sample_category_mapping(counts_df: pd.DataFrame, meta_df: pd.DataFrame, condition_col: str):
    print(f"Sample/category mapping used for DE (column='{condition_col}'):")
    if counts_df.empty:
        print("  (no samples available)")
        return
    if condition_col not in meta_df.columns:
        print(f"  (metadata column '{condition_col}' not found)")
        for sample in counts_df.index:
            print(f"  - sample={sample} category=NA")
        return
    for sample in counts_df.index:
        category = meta_df.loc[sample, condition_col] if sample in meta_df.index else "NA"
        print(f"  - sample={sample} category={category}")


def write_pca(
    counts_df: pd.DataFrame,
    meta_df: Optional[pd.DataFrame],
    condition_col: Optional[str],
    out_dir: Path,
    resume: bool,
):
    import matplotlib.pyplot as plt
    import numpy as np

    pca_tsv = out_dir / "pca_coordinates.tsv"
    pca_png = out_dir / "pca_plot.png"
    if should_skip(pca_tsv, resume, "PCA coordinates") and should_skip(pca_png, resume, "PCA plot"):
        return

    vst_like = np.log2(counts_df + 1.0)
    x = vst_like.to_numpy(dtype=float)
    x = x - x.mean(axis=0, keepdims=True)
    u, s, vt = np.linalg.svd(x, full_matrices=False)
    pc_scores = u[:, :2] * s[:2]
    var_exp = (s ** 2) / max((x.shape[0] - 1), 1)
    var_exp_ratio = var_exp / max(var_exp.sum(), 1e-12)
    use_groups = bool(meta_df is not None and condition_col and condition_col in meta_df.columns)
    pca_df = pd.DataFrame(
        {
            "sample": counts_df.index,
            "PC1": pc_scores[:, 0] if pc_scores.shape[1] > 0 else 0.0,
            "PC2": pc_scores[:, 1] if pc_scores.shape[1] > 1 else 0.0,
            (condition_col if use_groups else "group"): (
                meta_df[condition_col].astype(str).values if use_groups else "all_samples"
            ),
        }
    )
    pca_df.to_csv(pca_tsv, sep="\t", index=False)

    group_col = condition_col if use_groups else "group"
    fig, ax = plt.subplots(figsize=(6, 5))
    for cond, group in pca_df.groupby(group_col):
        ax.scatter(group["PC1"], group["PC2"], label=str(cond), s=70, alpha=0.85)
        for _, row in group.iterrows():
            ax.annotate(row["sample"], (row["PC1"], row["PC2"]), fontsize=8, alpha=0.85)
    ax.set_xlabel(f"PC1 ({var_exp_ratio[0] * 100:.1f}% var)" if len(var_exp_ratio) > 0 else "PC1")
    ax.set_ylabel(f"PC2 ({var_exp_ratio[1] * 100:.1f}% var)" if len(var_exp_ratio) > 1 else "PC2")
    ax.set_title("PCA of samples (log2 counts + 1)")
    if use_groups:
        ax.legend()
    fig.tight_layout()
    fig.savefig(pca_png, dpi=180)
    plt.close(fig)
    print(f"PCA outputs written: {pca_tsv}, {pca_png}")


def write_heatmap(counts_df: pd.DataFrame, out_dir: Path, resume: bool, top_n: int = 50):
    import matplotlib.pyplot as plt
    import numpy as np

    heatmap_png = out_dir / "heatmap_top_variable_genes.png"
    if should_skip(heatmap_png, resume, "heatmap"):
        return
    if counts_df.empty or counts_df.shape[0] < 2 or counts_df.shape[1] < 2:
        print("Skipping heatmap: need at least 2 samples and 2 genes.")
        return

    x = np.log2(counts_df + 1.0)
    top_genes = x.var(axis=0).sort_values(ascending=False).head(max(top_n, 2)).index
    z = x[top_genes]
    z = (z - z.mean(axis=0)) / z.std(axis=0).replace(0, 1)

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(z.T.to_numpy(), aspect="auto", interpolation="nearest", cmap="viridis")
    ax.set_title("Top variable genes (z-scored log2 counts)")
    ax.set_xlabel("Samples")
    ax.set_ylabel("Genes")
    ax.set_xticks(range(len(z.index)))
    ax.set_xticklabels(z.index, rotation=90, fontsize=7)
    fig.colorbar(im, ax=ax, label="z-score")
    fig.tight_layout()
    fig.savefig(heatmap_png, dpi=180)
    plt.close(fig)
    print(f"Heatmap written: {heatmap_png}")


def write_volcano(deseq_df: pd.DataFrame, out_dir: Path, resume: bool):
    import matplotlib.pyplot as plt
    import numpy as np

    volcano_png = out_dir / "volcano_plot.png"
    if should_skip(volcano_png, resume, "volcano plot"):
        return

    plot_df = deseq_df.copy()
    if "padj" not in plot_df.columns or "log2FoldChange" not in plot_df.columns:
        print("Skipping volcano plot: DESeq2 columns 'padj'/'log2FoldChange' not found.")
        return
    plot_df["padj"] = pd.to_numeric(plot_df["padj"], errors="coerce")
    plot_df["log2FoldChange"] = pd.to_numeric(plot_df["log2FoldChange"], errors="coerce")
    plot_df = plot_df.dropna(subset=["padj", "log2FoldChange"])
    if plot_df.empty:
        print("Skipping volcano plot: no valid DE rows.")
        return
    plot_df["neg_log10_padj"] = -np.log10(plot_df["padj"].clip(lower=1e-300))
    sig = (plot_df["padj"] < 0.05) & (plot_df["log2FoldChange"].abs() >= 1.0)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(plot_df.loc[~sig, "log2FoldChange"], plot_df.loc[~sig, "neg_log10_padj"], s=10, alpha=0.4, color="grey")
    ax.scatter(plot_df.loc[sig, "log2FoldChange"], plot_df.loc[sig, "neg_log10_padj"], s=14, alpha=0.7, color="crimson")
    ax.axhline(-math.log10(0.05), linestyle="--", color="black", linewidth=1)
    ax.axvline(-1.0, linestyle="--", color="black", linewidth=1)
    ax.axvline(1.0, linestyle="--", color="black", linewidth=1)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("Volcano plot")
    fig.tight_layout()
    fig.savefig(volcano_png, dpi=180)
    plt.close(fig)
    print(f"Volcano plot written: {volcano_png}")


def run_gsea(
    deseq_df: pd.DataFrame,
    gsea_gmt: str,
    out_dir: Path,
    resume: bool,
    gsea_min_size: int,
    gsea_max_size: int,
):
    gsea_dir = out_dir / "gsea"
    ensure_dir(gsea_dir)
    marker = gsea_dir / "gseapy.gene_set.prerank.report.csv"
    if should_skip(marker, resume, "GSEA report"):
        return
    try:
        import gseapy as gp
    except ImportError:
        print("Skipping GSEA: gseapy is not installed. Add it to environment to enable GSEA.")
        return

    rank_col = "stat" if "stat" in deseq_df.columns else "log2FoldChange"
    if rank_col not in deseq_df.columns:
        print("Skipping GSEA: no suitable ranking column ('stat' or 'log2FoldChange').")
        return
    ranking = deseq_df[["gene_id", rank_col]].copy()
    ranking[rank_col] = pd.to_numeric(ranking[rank_col], errors="coerce")
    ranking = ranking.dropna().sort_values(rank_col, ascending=False)
    if ranking.empty:
        print("Skipping GSEA: ranking table is empty.")
        return
    gp.prerank(
        rnk=ranking,
        gene_sets=gsea_gmt,
        outdir=str(gsea_dir),
        min_size=gsea_min_size,
        max_size=gsea_max_size,
        seed=42,
        verbose=True,
    )
    print(f"GSEA outputs written: {gsea_dir}")


def run_deseq(
    counts_matrix: Path,
    metadata_path: Path,
    condition_col: str,
    out_dir: Path,
    resume: bool,
    gsea_gmt: str,
    gsea_min_size: int,
    gsea_max_size: int,
):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    try:
        counts_df, meta_df, condition_col = _sample_ordered_counts_and_metadata(
            counts_matrix, metadata_path, condition_col
        )
    except Exception as exc:
        print(f"Skipping DE/volcano because metadata matching failed: {exc}")
        return

    # Ensure exact index compatibility for AnnData/PyDESeq2 internals.
    counts_df = counts_df.copy()
    meta_df = meta_df.copy()
    counts_df.index = pd.Index(counts_df.index.astype(str), name="sample")
    meta_df.index = pd.Index(meta_df.index.astype(str), name="sample")
    counts_df.columns = counts_df.columns.astype(str)
    counts_df = counts_df.loc[meta_df.index]
    meta_df = meta_df.loc[counts_df.index]
    print_sample_category_mapping(counts_df, meta_df, condition_col)

    out = out_dir / "deseq2_results.tsv"
    try:
        if should_skip(out, resume, "DESeq2 results"):
            res = pd.read_csv(out, sep="\t")
        else:
            dds = DeseqDataSet(counts=counts_df, metadata=meta_df, design=f"~{condition_col}")
            dds.deseq2()
            levels = [str(v) for v in pd.Series(meta_df[condition_col]).dropna().astype(str).unique().tolist()]
            if len(levels) < 2:
                raise ValueError(
                    f"Need at least 2 groups in metadata column '{condition_col}' for DESeq2; got {levels}."
                )
            contrast = [condition_col, levels[1], levels[0]]
            print(
                "DESeq2 contrast selected: "
                f"{condition_col}: {levels[1]} vs {levels[0]}"
            )
            stats = DeseqStats(dds, contrast=contrast, n_cpus=1)
            stats.summary()
            res = stats.results_df.reset_index().rename(columns={"index": "gene_id"})
            res.to_csv(out, sep="\t", index=False)
            print(f"Differential expression results written: {out}")
    except Exception as exc:
        print(f"Skipping DE/volcano because DESeq2 failed: {exc}")
        return

    write_pca(counts_df, meta_df, condition_col, out_dir, resume)
    write_heatmap(counts_df, out_dir, resume)
    try:
        write_volcano(res, out_dir, resume)
    except Exception as exc:
        print(f"Skipping volcano plot: {exc}")
    if gsea_gmt:
        try:
            run_gsea(res, gsea_gmt, out_dir, resume, gsea_min_size, gsea_max_size)
        except Exception as exc:
            print(f"Skipping GSEA: {exc}")


def build_parser():
    parser = argparse.ArgumentParser(description="QuickTranscriptome: minimal RNA-seq workflow")
    sub = parser.add_subparsers(dest="command", required=True)

    run_p = sub.add_parser("run", help="Run alignment + feature counting (+ optional DE)")
    run_p.add_argument("--reads-dir", required=True, help="Directory with *.fastq.gz files")
    run_p.add_argument("--species", default="parapsilosis", help="Species key (default: parapsilosis)")
    run_p.add_argument("--fasta-url", default=None, help="Custom FASTA URL")
    run_p.add_argument("--gff-url", default=None, help="Custom GFF3 URL")
    run_p.add_argument("--out-dir", default="results", help="Output directory")
    run_p.add_argument("--threads", type=int, default=4, help="Number of threads")
    run_p.add_argument(
        "--metadata",
        default=None,
        help="Optional sample metadata TSV/CSV with columns: sample,<condition>",
    )
    run_p.add_argument(
        "--condition-column",
        default=None,
        help="Condition/treatment column for DE (requires --metadata)",
    )
    run_p.add_argument(
        "--resume",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Resume from existing outputs if present (default: True). Use --no-resume to force rerun.",
    )
    run_p.add_argument(
        "--gsea-gmt",
        default=None,
        help="Optional path to GMT gene-set file for preranked GSEA.",
    )
    run_p.add_argument("--gsea-min-size", type=int, default=10, help="Minimum gene-set size for GSEA.")
    run_p.add_argument("--gsea-max-size", type=int, default=500, help="Maximum gene-set size for GSEA.")
    run_p.set_defaults(func=align_and_count)
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
