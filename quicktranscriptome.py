#!/usr/bin/env python3
import argparse
import base64
import gzip
import html
import json
import math
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlretrieve, urlopen

import pandas as pd


def normalize_gene_identifier(raw: str) -> str:
    """Match GFF/featureCounts IDs to CGD GAF symbols (strip Ensembl-style gene: prefix)."""
    s = str(raw).strip()
    if s.lower().startswith("gene:"):
        return s[5:].lstrip()
    return s


def _fetch_go_term_names(go_ids: list[str]) -> dict[str, str]:
    """Map GO:xxxxxx to term names via QuickGO (best effort; offline returns {})."""
    uniq = sorted({g.strip() for g in go_ids if str(g).strip().startswith("GO:")})
    if not uniq:
        return {}
    out: dict[str, str] = {}
    chunk = 40
    for i in range(0, len(uniq), chunk):
        part = uniq[i : i + chunk]
        url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/" + ",".join(part)
        try:
            req = Request(url, headers={"Accept": "application/json"})
            with urlopen(req, timeout=45) as resp:
                data = json.loads(resp.read().decode())
            for t in data.get("results", []):
                tid = t.get("id")
                name = t.get("name")
                if tid and name:
                    out[str(tid)] = str(name)
        except (HTTPError, URLError, OSError, TimeoutError, ValueError, json.JSONDecodeError):
            pass
    return out


def _format_go_axis_label(term: str, name_by_id: dict[str, str]) -> str:
    t = str(term).strip()
    if t.startswith("GO:") and t in name_by_id:
        return f"{name_by_id[t]} ({t})"
    return t


def _truncate_label(s: str, max_len: int = 72) -> str:
    s = str(s)
    return s if len(s) <= max_len else s[: max_len - 1] + "…"


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
        "gsea_gaf_urls": [
            "http://www.candidagenome.org/download/go/cgd_C_parapsilosis_CDC317.gaf.gz",
            "http://www.candidagenome.org/download/go/gene_association.cgd.gz",
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


def _build_gmt_from_gaf(gaf_path: Path, gmt_path: Path) -> Path:
    """Convert a GAF file into a simple GO-term GMT file."""
    go_to_genes = defaultdict(set)
    open_func = gzip.open if gaf_path.suffix == ".gz" else open
    with open_func(gaf_path, "rt", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.strip() or line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            qualifier = parts[3]
            if "NOT" in qualifier.split("|"):
                continue
            go_id = parts[4].strip()
            if not go_id:
                continue
            # CGD: col2 = CAL object ID, col3 = systematic gene symbol (CPAR2_*).
            # Count matrices use the same symbols as GFF gene IDs, not CAL IDs.
            symbol = parts[2].strip() if len(parts) > 2 else ""
            obj_id = parts[1].strip() if len(parts) > 1 else ""
            gene_id = normalize_gene_identifier(symbol or obj_id)
            if not gene_id:
                continue
            go_to_genes[go_id].add(gene_id)

    with open(gmt_path, "w", encoding="utf-8") as out:
        for go_id in sorted(go_to_genes):
            genes = sorted(go_to_genes[go_id])
            if not genes:
                continue
            out.write(f"{go_id}\tCGD_GO\t" + "\t".join(genes) + "\n")
    return gmt_path


def resolve_gsea_gmt(
    species: str,
    refs_dir: Path,
    user_gsea_gmt: Optional[str],
    user_gsea_gaf_url: Optional[str],
    resume: bool,
) -> Optional[str]:
    if user_gsea_gmt:
        return user_gsea_gmt

    species_defaults = SPECIES_DEFAULTS.get(species, {})
    gaf_urls = normalize_urls(user_gsea_gaf_url) or normalize_urls(species_defaults.get("gsea_gaf_urls"))
    if not gaf_urls:
        return None

    gsea_ref_dir = refs_dir / "gsea"
    ensure_dir(gsea_ref_dir)
    # Suffix distinguishes symbol-based GMT from older CAL-ID caches.
    gmt_out = gsea_ref_dir / f"{species}_go_from_gaf_symbol.gmt"
    if should_skip(gmt_out, resume, "auto-generated GSEA GMT"):
        return str(gmt_out)

    gaf_path = download_first_available(gaf_urls, gsea_ref_dir)
    _build_gmt_from_gaf(gaf_path, gmt_out)
    print(f"Auto-generated GSEA GMT from GAF: {gmt_out}")
    return str(gmt_out)


def _summary_plot_candidates(counts_dir: Path) -> list[tuple[str, Path, str]]:
    """Title, path, short description for report."""
    return [
        (
            "PCA",
            counts_dir / "pca_plot.png",
            "Principal components on log2(counts+1), mean-centered. Samples separate in PC space when "
            "replicates of the same condition cluster together.",
        ),
        (
            "Top variable genes",
            counts_dir / "heatmap_top_variable_genes.png",
            "Genes with the highest variance across samples (after log2 transform), z-scored per gene. "
            "Highlights global sample-to-sample heterogeneity before differential testing.",
        ),
        (
            "Volcano (DESeq2)",
            counts_dir / "volcano_plot.png",
            "Effect size (log2 fold change) versus significance (-log10 adjusted p-value). "
            "Labels mark the 20 genes with smallest adjusted p-value. "
            "Horizontal lines mark typical padj and |log2FC| cutoffs used for visualization.",
        ),
        (
            "Top DE genes",
            counts_dir / "heatmap_top_de_genes.png",
            "Top up- and down-regulated genes from DESeq2 (see contrast in title), z-scored log2 counts. "
            "Rows: genes; columns: samples ordered by condition.",
        ),
        (
            "GSEA prerank (NES)",
            counts_dir / "gsea" / "gsea_prerank_nes_dotplot.png",
            "Gene-set enrichment using ranks from DESeq2. Shows the five most negative and five most "
            "positive NES terms; y-axis labels use GO names from QuickGO when online. "
            "Marker size reflects significance (-log10 p). Red = positive NES, blue = negative.",
        ),
    ]


def _collect_run_summary_stats(out_dir: Path, counts_dir: Path, args: argparse.Namespace) -> list[tuple[str, str]]:
    """Key/value rows for the summary table."""
    rows: list[tuple[str, str]] = []
    cm = counts_dir / "counts_matrix.tsv"
    if cm.exists():
        try:
            header = cm.read_text(encoding="utf-8", errors="replace").splitlines()[:1]
            if header:
                ncols = header[0].count("\t") + 1
                n_samples = max(0, ncols - 1)
                rows.append(("Samples in matrix", str(n_samples)))
            n_lines = sum(1 for _ in cm.open("r", encoding="utf-8", errors="replace")) - 1
            rows.append(("Genes in count matrix", str(max(0, n_lines))))
        except OSError:
            pass

    de_path = counts_dir / "deseq2_results.tsv"
    if de_path.exists():
        try:
            de = pd.read_csv(de_path, sep="\t")
            if "padj" in de.columns and "log2FoldChange" in de.columns:
                padj = pd.to_numeric(de["padj"], errors="coerce")
                lfc = pd.to_numeric(de["log2FoldChange"], errors="coerce")
                rows.append(("DE genes (padj < 0.05)", str(int(((padj < 0.05) & lfc.notna()).sum()))))
                rows.append(("DE genes (padj < 0.01)", str(int(((padj < 0.01) & lfc.notna()).sum()))))
        except (OSError, ValueError):
            pass

    gsea_rep = counts_dir / "gsea" / "gseapy.gene_set.prerank.report.csv"
    if gsea_rep.exists():
        try:
            gr = pd.read_csv(gsea_rep, sep=",")
            rows.append(("GSEA gene sets in report", str(len(gr))))
            if "NES" in gr.columns:
                nes = pd.to_numeric(gr["NES"], errors="coerce")
                rows.append(("GSEA sets |NES| max", f"{nes.abs().max():.3f}" if nes.notna().any() else "—"))
        except (OSError, ValueError):
            pass

    rows.append(("BAM directory", str((out_dir / "bam").resolve())))
    rows.append(("Reference index", str((out_dir / "refs" / "reference.mmi").resolve())))
    return rows


def write_run_summary(out_dir: Path, counts_dir: Path, args: argparse.Namespace, resume: bool):
    """Write multi-section HTML and multi-page PDF with figures, table, and short descriptions."""
    html_path = out_dir / "run_summary.html"
    pdf_path = out_dir / "run_summary.pdf"
    if resume and html_path.exists() and pdf_path.exists():
        print(f"Resume enabled: using existing run summary: {html_path}, {pdf_path}")
        return

    plot_sections = [(t, p, d) for t, p, d in _summary_plot_candidates(counts_dir) if p.is_file()]
    stats_rows = _collect_run_summary_stats(out_dir, counts_dir, args)

    def _img_data_uri(path: Path) -> str:
        raw = path.read_bytes()
        b64 = base64.standard_b64encode(raw).decode("ascii")
        return f"data:image/png;base64,{b64}"

    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    lines = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="utf-8">',
        "<title>QuickTranscriptome run summary</title>",
        "<style>",
        "body { font-family: system-ui, Segoe UI, sans-serif; margin: 16px 20px; font-size: 14px; color: #1a1a1a; line-height: 1.5; max-width: 960px; }",
        "h1 { font-size: 1.35rem; margin: 0 0 12px 0; }",
        "h2 { font-size: 1.05rem; margin: 28px 0 8px 0; color: #111; border-bottom: 1px solid #ddd; padding-bottom: 4px; }",
        "p.desc { margin: 8px 0 12px 0; color: #333; }",
        ".meta { background: #f4f4f5; padding: 12px 14px; border-radius: 8px; margin-bottom: 20px; line-height: 1.5; }",
        ".meta code { font-size: 12px; word-break: break-all; }",
        "table.stats { border-collapse: collapse; width: 100%; max-width: 640px; margin: 12px 0 24px 0; font-size: 13px; }",
        "table.stats th, table.stats td { border: 1px solid #ccc; padding: 8px 10px; text-align: left; }",
        "table.stats th { background: #eee; font-weight: 600; }",
        "table.stats tr:nth-child(even) { background: #fafafa; }",
        ".fig { border: 1px solid #ddd; border-radius: 6px; padding: 10px; background: #fff; margin-bottom: 8px; }",
        ".fig img { width: 100%; height: auto; display: block; object-fit: contain; }",
        ".fig:not(.gsea) img { max-height: 420px; }",
        ".fig.gsea img { max-height: none; }",
        "@media print { .section { page-break-after: always; } .section:last-child { page-break-after: auto; } }",
        "</style>",
        "</head>",
        "<body>",
        "<h1>QuickTranscriptome — run summary</h1>",
        '<section class="section meta">',
        f"<p><strong>Generated</strong>: {html.escape(now)}</p>",
        f"<p><strong>Species</strong>: {html.escape(str(args.species))} &nbsp;|&nbsp; "
        f"<strong>Threads</strong>: {args.threads}</p>",
        f"<p><strong>Reads</strong><br><code>{html.escape(str(Path(args.reads_dir).resolve()))}</code></p>",
        f"<p><strong>Output</strong><br><code>{html.escape(str(out_dir.resolve()))}</code></p>",
    ]
    if args.metadata:
        lines.append(
            f"<p><strong>Metadata</strong><br><code>{html.escape(str(Path(args.metadata).resolve()))}</code><br>"
            f"Condition column: {html.escape(str(args.condition_column))}</p>"
        )
    else:
        lines.append("<p><strong>Differential expression</strong>: not run (no metadata).</p>")
    lines.append(
        f"<p><strong>GSEA gene sets (GMT)</strong>: {html.escape(str(args.gsea_gmt or 'auto from species defaults if available'))}</p>"
    )
    lines.append("</section>")

    lines.append('<section class="section">')
    lines.append("<h2>Run statistics</h2>")
    lines.append("<p class=\"desc\">Summary metrics from this output folder (counts, DESeq2, GSEA when present).</p>")
    if stats_rows:
        lines.append('<table class="stats">')
        lines.append("<tr><th>Metric</th><th>Value</th></tr>")
        for k, v in stats_rows:
            lines.append(f"<tr><td>{html.escape(k)}</td><td>{html.escape(str(v))}</td></tr>")
        lines.append("</table>")
    else:
        lines.append("<p>No statistics could be loaded (run may still be in progress).</p>")
    lines.append("</section>")

    if not plot_sections:
        lines.append('<section class="section"><h2>Figures</h2><p>No plot PNGs found under <code>counts/</code> yet.</p></section>')
    else:
        for title, path, blurb in plot_sections:
            cls = "fig gsea" if "GSEA" in title else "fig"
            lines.append('<section class="section">')
            lines.append(f"<h2>{html.escape(title)}</h2>")
            lines.append(f'<p class="desc">{html.escape(blurb)}</p>')
            lines.append(f'<div class="{cls}">')
            lines.append(f'<img src="{_img_data_uri(path)}" alt="{html.escape(title)}">')
            lines.append("</div></section>")

    lines.extend(["</body>", "</html>"])
    html_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Run summary HTML written: {html_path}")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.image as mpimg
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError as exc:
        print(f"Skipping run summary PDF (matplotlib unavailable): {exc}")
        return

    def _wrap_text(s: str, width: int = 95) -> str:
        words = s.split()
        lines_out: list[str] = []
        cur: list[str] = []
        n = 0
        for w in words:
            if n + len(w) + len(cur) > width and cur:
                lines_out.append(" ".join(cur))
                cur = [w]
                n = len(w)
            else:
                cur.append(w)
                n = sum(len(x) + 1 for x in cur)
        if cur:
            lines_out.append(" ".join(cur))
        return "\n".join(lines_out)

    with PdfPages(pdf_path) as pdf:
        fig = plt.figure(figsize=(8.5, 11))
        fig.patch.set_facecolor("white")
        y = 0.97
        fig.text(0.08, y, "QuickTranscriptome — run summary", fontsize=15, fontweight="bold")
        y -= 0.045
        fig.text(0.08, y, f"Generated: {now}", fontsize=9)
        y -= 0.028
        fig.text(0.08, y, f"Species: {args.species}  |  Threads: {args.threads}", fontsize=9)
        y -= 0.035
        fig.text(0.08, y, "Reads:", fontsize=9, fontweight="bold")
        y -= 0.022
        fig.text(0.08, y, str(Path(args.reads_dir).resolve()), fontsize=7, family="monospace")
        y -= 0.04
        fig.text(0.08, y, "Output:", fontsize=9, fontweight="bold")
        y -= 0.022
        fig.text(0.08, y, str(out_dir.resolve()), fontsize=7, family="monospace")
        y -= 0.045
        if stats_rows:
            tbl_h = min(0.52, 0.08 + 0.034 * len(stats_rows))
            ax_tbl = fig.add_axes([0.08, 0.28, 0.84, tbl_h])
            ax_tbl.axis("off")
            tbl = ax_tbl.table(
                cellText=[[a, b] for a, b in stats_rows],
                colLabels=["Metric", "Value"],
                loc="upper center",
                cellLoc="left",
            )
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(8)
            tbl.scale(1.0, 1.4)
        else:
            fig.text(0.08, 0.45, "No statistics table (outputs missing or incomplete).", fontsize=10)

        pdf.savefig(fig)
        plt.close(fig)

        if not plot_sections:
            fig2 = plt.figure(figsize=(8.5, 11))
            fig2.text(0.5, 0.5, "No plot PNGs were found for this run.", ha="center", va="center", fontsize=12)
            pdf.savefig(fig2)
            plt.close(fig2)
        else:
            for title, path, blurb in plot_sections:
                figp = plt.figure(figsize=(8.5, 11))
                figp.patch.set_facecolor("white")
                figp.text(0.07, 0.98, title, fontsize=13, fontweight="bold", va="top")
                blurb_short = "\n".join(_wrap_text(blurb, 92).split("\n")[:8])
                figp.text(
                    0.07,
                    0.91,
                    blurb_short,
                    fontsize=8.5,
                    va="top",
                    ha="left",
                    linespacing=1.22,
                )
                img = mpimg.imread(str(path))
                ax_img = figp.add_axes([0.06, 0.06, 0.88, 0.72])
                ax_img.imshow(img, aspect="auto", interpolation="bilinear")
                ax_img.axis("off")
                pdf.savefig(figp)
                plt.close(figp)

    print(f"Run summary PDF written: {pdf_path}")


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
    clean["gene_id"] = clean["gene_id"].map(normalize_gene_identifier)
    clean_out = counts_dir / "counts_matrix.tsv"
    if not should_skip(clean_out, args.resume, "counts matrix"):
        clean.to_csv(clean_out, sep="\t", index=False)
        print(f"Counts matrix written: {clean_out}")

    counts_for_plots = clean.set_index("gene_id").T.astype(int)
    gsea_gmt = resolve_gsea_gmt(
        species=args.species,
        refs_dir=refs_dir,
        user_gsea_gmt=args.gsea_gmt,
        user_gsea_gaf_url=args.gsea_gaf_url,
        resume=args.resume,
    )
    if args.metadata and args.condition_column:
        run_deseq(
            clean_out,
            Path(args.metadata),
            args.condition_column,
            counts_dir,
            args.resume,
            gsea_gmt,
            args.gsea_min_size,
            args.gsea_max_size,
            args.gsea_dotplot_max_terms,
            args.de_heatmap_top_per,
            args.de_heatmap_padj,
        )
    else:
        write_pca(counts_for_plots, None, None, counts_dir, args.resume)
        write_heatmap(counts_for_plots, counts_dir, args.resume)

    write_run_summary(out_dir, counts_dir, args, args.resume)


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


def write_heatmap(
    counts_df: pd.DataFrame,
    out_dir: Path,
    resume: bool,
    top_n: int = 50,
    meta_df: Optional[pd.DataFrame] = None,
    condition_col: Optional[str] = None,
):
    import matplotlib.pyplot as plt
    import numpy as np

    heatmap_png = out_dir / "heatmap_top_variable_genes.png"
    if should_skip(heatmap_png, resume, "heatmap"):
        return
    if counts_df.empty or counts_df.shape[0] < 2 or counts_df.shape[1] < 2:
        print("Skipping heatmap: need at least 2 samples and 2 genes.")
        return

    x = np.log2(counts_df + 1.0)
    if meta_df is not None and condition_col and condition_col in meta_df.columns:
        common = x.index.intersection(meta_df.index)
        if len(common) >= 2:
            order = (
                meta_df.loc[common, condition_col]
                .astype(str)
                .sort_values()
                .index
            )
            x = x.loc[order]

    top_genes = x.var(axis=0).sort_values(ascending=False).head(max(top_n, 2)).index
    z = x[top_genes]
    z = (z - z.mean(axis=0)) / z.std(axis=0).replace(0, 1)

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(z.T.to_numpy(), aspect="auto", interpolation="nearest", cmap="viridis")
    ax.set_title("Top variable genes (z-scored log2 counts)")
    ax.set_xlabel("Samples")
    ax.set_ylabel("Genes")
    ax.set_xticks(range(len(z.index)))
    if meta_df is not None and condition_col and condition_col in meta_df.columns:
        sample_labels = [f"{s} ({meta_df.loc[s, condition_col]})" for s in z.index]
        ax.set_xticklabels(sample_labels, rotation=90, fontsize=7)
    else:
        ax.set_xticklabels(z.index, rotation=90, fontsize=7)
    fig.colorbar(im, ax=ax, label="z-score")
    fig.tight_layout()
    fig.savefig(heatmap_png, dpi=180)
    plt.close(fig)
    print(f"Heatmap written: {heatmap_png}")


def write_de_top_heatmap(
    counts_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    condition_col: str,
    deseq_df: pd.DataFrame,
    out_dir: Path,
    resume: bool,
    top_per_direction: int,
    padj_cutoff: float,
    numerator_level: str,
    denominator_level: str,
):
    """Heatmap of log2 counts (z-scored per gene) for top genes up/down in the DESeq2 contrast."""
    import matplotlib.pyplot as plt
    import numpy as np

    heatmap_png = out_dir / "heatmap_top_de_genes.png"
    if should_skip(heatmap_png, resume, "DE top genes heatmap"):
        return
    if counts_df.empty or counts_df.shape[0] < 2:
        print("Skipping DE heatmap: need at least 2 samples.")
        return
    if "padj" not in deseq_df.columns or "log2FoldChange" not in deseq_df.columns:
        print("Skipping DE heatmap: DESeq2 columns 'padj'/'log2FoldChange' not found.")
        return

    de = deseq_df.copy()
    de["gene_id"] = de["gene_id"].astype(str).str.strip().map(normalize_gene_identifier)
    de["padj"] = pd.to_numeric(de["padj"], errors="coerce")
    de["log2FoldChange"] = pd.to_numeric(de["log2FoldChange"], errors="coerce")
    de = de.dropna(subset=["log2FoldChange"])

    # Map normalized gene id -> actual column name in counts (PyDESeq2/CSV may differ slightly from matrix cols).
    col_by_norm: dict[str, str] = {}
    for c in counts_df.columns:
        key = normalize_gene_identifier(str(c))
        if key not in col_by_norm:
            col_by_norm[key] = str(c)

    def _to_count_col(g) -> Optional[str]:
        if g is None or pd.isna(g):
            return None
        key = normalize_gene_identifier(str(g).strip())
        if key in col_by_norm:
            return col_by_norm[key]
        gs = str(g).strip()
        if gs in counts_df.columns:
            return gs
        return None

    sig = de["padj"].notna() & (de["padj"] < padj_cutoff)
    pos = de.loc[sig & (de["log2FoldChange"] > 0)].sort_values("log2FoldChange", ascending=False)
    neg = de.loc[sig & (de["log2FoldChange"] < 0)].sort_values("log2FoldChange", ascending=True)

    if len(pos) < top_per_direction:
        pos = de.loc[de["log2FoldChange"] > 0].sort_values("log2FoldChange", ascending=False).head(top_per_direction)
    else:
        pos = pos.head(top_per_direction)

    if len(neg) < top_per_direction:
        neg = de.loc[de["log2FoldChange"] < 0].sort_values("log2FoldChange", ascending=True).head(top_per_direction)
    else:
        neg = neg.head(top_per_direction)

    gene_ids_raw = list(neg["gene_id"]) + list(pos["gene_id"])
    gene_ids: list[str] = []
    seen: set[str] = set()
    for g in gene_ids_raw:
        col = _to_count_col(g)
        if col is None or col in seen:
            continue
        seen.add(col)
        gene_ids.append(col)

    if len(gene_ids) < 2:
        pool = de.copy()
        pool["_abs_lfc"] = pool["log2FoldChange"].abs()
        pool = pool.sort_values("_abs_lfc", ascending=False)
        for _, row in pool.iterrows():
            col = _to_count_col(row["gene_id"])
            if col is None or col in seen:
                continue
            seen.add(col)
            gene_ids.append(col)
            if len(gene_ids) >= max(2, min(2 * top_per_direction, 50)):
                break

    if len(gene_ids) < 2:
        ex_de = de["gene_id"].head(3).tolist()
        ex_ct = [normalize_gene_identifier(str(c)) for c in counts_df.columns[:3]]
        print(
            "Skipping DE heatmap: fewer than two genes from DESeq2 results could be matched to count-matrix "
            f"columns (check gene_id consistency). Example DE gene_ids: {ex_de}; "
            f"example count columns (normalized): {ex_ct}."
        )
        return

    x = np.log2(counts_df + 1.0)
    if meta_df is not None and condition_col and condition_col in meta_df.columns:
        common = x.index.intersection(meta_df.index)
        if len(common) >= 2:
            order = (
                meta_df.loc[common, condition_col]
                .astype(str)
                .sort_values()
                .index
            )
            x = x.loc[order]

    z = x[gene_ids]
    z = (z - z.mean(axis=0)) / z.std(axis=0).replace(0, 1)

    fig_h = min(max(5.0, 0.22 * len(gene_ids) + 2.0), 100.0)
    fig, ax = plt.subplots(figsize=(10, fig_h))
    arr = z.T.to_numpy()
    vmax = float(np.nanpercentile(np.abs(arr), 99)) if np.isfinite(arr).any() else 3.0
    vmax = max(vmax, 0.5)
    im = ax.imshow(
        arr,
        aspect="auto",
        interpolation="nearest",
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_title(
        f"Top DE genes (z-scored log2 counts)\n"
        f"positive log2FC = higher in {numerator_level} vs {denominator_level}"
    )
    ax.set_xlabel("Samples")
    ax.set_ylabel("Genes (downregulated, then upregulated)")
    ax.set_xticks(range(len(z.index)))
    if meta_df is not None and condition_col and condition_col in meta_df.columns:
        sample_labels = [f"{s} ({meta_df.loc[s, condition_col]})" for s in z.index]
        ax.set_xticklabels(sample_labels, rotation=90, fontsize=7)
    else:
        ax.set_xticklabels(z.index, rotation=90, fontsize=7)
    ax.set_yticks(range(len(gene_ids)))
    ax.set_yticklabels(gene_ids, fontsize=max(5, 8 - len(gene_ids) // 40))
    fig.colorbar(im, ax=ax, label="z-score (per gene)")
    fig.tight_layout()
    fig.savefig(heatmap_png, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"DE top genes heatmap written: {heatmap_png}")


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

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(plot_df.loc[~sig, "log2FoldChange"], plot_df.loc[~sig, "neg_log10_padj"], s=10, alpha=0.4, color="grey")
    ax.scatter(plot_df.loc[sig, "log2FoldChange"], plot_df.loc[sig, "neg_log10_padj"], s=14, alpha=0.7, color="crimson")
    ax.axhline(-math.log10(0.05), linestyle="--", color="black", linewidth=1)
    ax.axvline(-1.0, linestyle="--", color="black", linewidth=1)
    ax.axvline(1.0, linestyle="--", color="black", linewidth=1)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title("Volcano plot (labels: top 20 genes by smallest adjusted p-value)")
    label_df = plot_df.sort_values("padj", ascending=True).head(20)
    for _, row in label_df.iterrows():
        gid = str(row["gene_id"])[:40]
        ax.annotate(
            gid,
            (row["log2FoldChange"], row["neg_log10_padj"]),
            fontsize=5.5,
            alpha=0.92,
            ha="left",
            va="bottom",
            bbox=dict(boxstyle="round,pad=0.12", facecolor="white", edgecolor="0.7", linewidth=0.3, alpha=0.85),
        )
    fig.tight_layout()
    fig.savefig(volcano_png, dpi=180)
    plt.close(fig)
    print(f"Volcano plot written: {volcano_png}")


def write_gsea_nes_dotplot(
    report_path: Path,
    out_path: Path,
    resume: bool,
    max_terms: int,
):
    """Five lowest- and five highest-NES terms; marker area scales with significance (-log10 p)."""
    if max_terms:
        print(
            "Note: --gsea-dotplot-max-terms is ignored; the NES dot plot shows "
            "5 lowest + 5 highest NES terms."
        )
    if should_skip(out_path, resume, "GSEA NES dot plot"):
        return
    if not report_path.exists():
        print("Skipping GSEA dot plot: report CSV not found.")
        return

    import matplotlib.pyplot as plt
    import numpy as np

    df = pd.read_csv(report_path)
    if df.empty:
        print("Skipping GSEA dot plot: empty report.")
        return

    term_col = "Term" if "Term" in df.columns else df.columns[0]
    if "NES" not in df.columns:
        print("Skipping GSEA dot plot: no NES column in report.")
        return

    nes = pd.to_numeric(df["NES"], errors="coerce")
    terms = df[term_col].astype(str)
    df_plot = pd.DataFrame({"term": terms, "nes": nes}).dropna(subset=["nes"])
    if df_plot.empty:
        print("Skipping GSEA dot plot: no valid NES values.")
        return

    df_plot = df_plot.drop_duplicates(subset=["term"], keep="first")
    neg_pool = df_plot[df_plot["nes"] < 0]
    pos_pool = df_plot[df_plot["nes"] > 0]
    neg5 = neg_pool.nsmallest(min(5, len(neg_pool)), "nes") if len(neg_pool) else pd.DataFrame()
    pos5 = pos_pool.nlargest(min(5, len(pos_pool)), "nes") if len(pos_pool) else pd.DataFrame()
    df_plot = pd.concat([neg5, pos5], ignore_index=True)
    if df_plot.empty:
        print("Skipping GSEA dot plot: no terms after selecting top/bottom NES.")
        return
    df_plot = df_plot.sort_values("nes", ascending=True).reset_index(drop=True)

    go_ids_for_lookup = df_plot["term"].astype(str).tolist()
    name_by_id = _fetch_go_term_names(go_ids_for_lookup)
    if any(str(t).startswith("GO:") for t in go_ids_for_lookup) and not name_by_id:
        print(
            "GSEA dot plot: could not fetch GO term names from QuickGO (offline or API error); "
            "using term IDs only."
        )

    p_col = next(
        (c for c in ("NOM p-val", "FDR q-val", "FWER p-val") if c in df.columns),
        None,
    )
    if p_col is not None:
        _pv = df[[term_col, p_col]].drop_duplicates(subset=[term_col])
        _pv[term_col] = _pv[term_col].astype(str)
        pmap = _pv.set_index(term_col)[p_col]
        pvals = df_plot["term"].map(pmap)
        pvals = pd.to_numeric(pvals, errors="coerce").fillna(1.0).clip(lower=1e-300)
        size_metric = -np.log10(pvals.to_numpy(dtype=float))
        size_label = f"-log10({p_col})"
    else:
        size_metric = df_plot["nes"].abs().to_numpy(dtype=float)
        size_label = "|NES| (no p-values in report)"
        print("GSEA dot plot: using |NES| for marker size (no p-value column found).")

    smin, smax = 18.0, 220.0
    sm = np.asarray(size_metric, dtype=float)
    valid = np.isfinite(sm)
    if not valid.any():
        sm = np.ones(len(df_plot)) * ((smin + smax) / 2)
    else:
        lo, hi = np.nanmin(sm[valid]), np.nanmax(sm[valid])
        if hi <= lo:
            sm = np.full(len(df_plot), (smin + smax) / 2)
        else:
            sm = np.where(valid, smin + (sm - lo) / (hi - lo) * (smax - smin), (smin + smax) / 2)

    y = np.arange(len(df_plot))
    n_t = len(df_plot)
    fig_w = 10.0
    fig_h = min(max(5.5, 0.19 * n_t + 2.0), 120.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    colors = np.where(df_plot["nes"].to_numpy() >= 0, "#c0392b", "#2980b9")
    ax.grid(True, axis="x", linestyle=":", alpha=0.45, zorder=0)
    ax.set_axisbelow(True)
    ax.scatter(
        df_plot["nes"],
        y,
        s=sm,
        c=colors,
        alpha=0.88,
        edgecolors="white",
        linewidths=0.55,
        zorder=2,
    )
    ax.axvline(0.0, color="0.35", linewidth=1.0, linestyle="--", zorder=1)
    ax.set_yticks(y)
    raw_terms = df_plot["term"].astype(str).tolist()
    labels = [_truncate_label(_format_go_axis_label(t, name_by_id)) for t in raw_terms]
    ytick_fs = max(7.0, min(9.5, 12.0 - n_t / 12.0))
    ax.set_yticklabels(labels, fontsize=ytick_fs)
    ax.set_xlabel("NES (normalized enrichment score)", fontsize=11)
    ax.set_ylabel("Gene set (GO term)", fontsize=11)
    ax.tick_params(axis="x", labelsize=10)
    ax.set_title(
        f"GSEA prerank (5 lowest + 5 highest NES; marker size ∝ {size_label})",
        fontsize=12,
        pad=10,
    )
    max_lab = max((len(str(l)) for l in labels), default=12)
    left_m = min(0.46, 0.11 + min(max_lab, 42) * 0.0085)
    fig.subplots_adjust(left=left_m, right=0.98, top=0.96, bottom=0.05)
    fig.savefig(out_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"GSEA NES dot plot written: {out_path}")


def run_gsea(
    deseq_df: pd.DataFrame,
    gsea_gmt: str,
    out_dir: Path,
    resume: bool,
    gsea_min_size: int,
    gsea_max_size: int,
    gsea_dotplot_max_terms: int,
):
    gsea_dir = out_dir / "gsea"
    ensure_dir(gsea_dir)
    report = gsea_dir / "gseapy.gene_set.prerank.report.csv"
    dotplot_png = gsea_dir / "gsea_prerank_nes_dotplot.png"

    skip_prerank = resume and report.exists()
    if skip_prerank:
        print(f"Resume enabled: using existing GSEA report: {report}")

    if not skip_prerank:
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
        ranking["gene_id"] = ranking["gene_id"].map(normalize_gene_identifier)
        ranking = ranking.drop_duplicates(subset=["gene_id"], keep="first")
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

    if not report.exists():
        print("Skipping GSEA dot plot: no GSEA report CSV.")
        return

    write_gsea_nes_dotplot(report, dotplot_png, resume, gsea_dotplot_max_terms)


def run_deseq(
    counts_matrix: Path,
    metadata_path: Path,
    condition_col: str,
    out_dir: Path,
    resume: bool,
    gsea_gmt: str,
    gsea_min_size: int,
    gsea_max_size: int,
    gsea_dotplot_max_terms: int,
    de_heatmap_top_per: int,
    de_heatmap_padj: float,
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
    counts_df = counts_df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    meta_df[condition_col] = meta_df[condition_col].astype(str).str.strip()
    meta_df.loc[meta_df[condition_col].isin(["", "nan", "None", "<NA>"]), condition_col] = pd.NA
    valid = meta_df[condition_col].notna()
    counts_df = counts_df.loc[valid]
    meta_df = meta_df.loc[valid]
    if counts_df.empty or len(meta_df) < 2:
        raise ValueError(
            f"Not enough valid samples after cleaning metadata column '{condition_col}' for DESeq2."
        )
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
    write_heatmap(counts_df, out_dir, resume, meta_df=meta_df, condition_col=condition_col)
    try:
        write_volcano(res, out_dir, resume)
    except Exception as exc:
        print(f"Skipping volcano plot: {exc}")
    levels_plot = [str(v) for v in pd.Series(meta_df[condition_col]).dropna().astype(str).unique().tolist()]
    num_l = levels_plot[1] if len(levels_plot) >= 2 else "?"
    den_l = levels_plot[0] if len(levels_plot) >= 2 else "?"
    try:
        write_de_top_heatmap(
            counts_df,
            meta_df,
            condition_col,
            res,
            out_dir,
            resume,
            de_heatmap_top_per,
            de_heatmap_padj,
            num_l,
            den_l,
        )
    except Exception as exc:
        print(f"Skipping DE top genes heatmap: {exc}")
    if gsea_gmt:
        try:
            run_gsea(
                res,
                gsea_gmt,
                out_dir,
                resume,
                gsea_min_size,
                gsea_max_size,
                gsea_dotplot_max_terms,
            )
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
        help="Optional path to GMT gene-set file for preranked GSEA. If omitted, an organism GAF may be auto-downloaded and converted to GMT when available.",
    )
    run_p.add_argument(
        "--gsea-gaf-url",
        default=None,
        help="Optional URL to a GAF file to auto-build GMT for GSEA (used when --gsea-gmt is not provided).",
    )
    run_p.add_argument("--gsea-min-size", type=int, default=10, help="Minimum gene-set size for GSEA.")
    run_p.add_argument("--gsea-max-size", type=int, default=500, help="Maximum gene-set size for GSEA.")
    run_p.add_argument(
        "--gsea-dotplot-max-terms",
        type=int,
        default=0,
        help="Ignored: NES dot plot is fixed at five lowest + five highest NES (kept for CLI compatibility).",
    )
    run_p.add_argument(
        "--de-heatmap-top-per",
        type=int,
        default=25,
        help="How many top up- and top down-regulated genes to show on the DE heatmap (each direction).",
    )
    run_p.add_argument(
        "--de-heatmap-padj",
        type=float,
        default=0.05,
        help="Adjusted p-value cutoff when choosing top DE genes for the DE heatmap.",
    )
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
