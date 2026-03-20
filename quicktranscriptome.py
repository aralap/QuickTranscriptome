#!/usr/bin/env python3
import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path
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
        if bam.exists():
            print(f"Using existing BAM: {bam}")
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
    run_cmd(
        [
            "featureCounts",
            "-T",
            str(args.threads),
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
    clean_out = counts_dir / "counts_matrix.tsv"
    clean.to_csv(clean_out, sep="\t", index=False)
    print(f"Counts matrix written: {clean_out}")

    if args.metadata and args.condition_column:
        run_deseq(clean_out, Path(args.metadata), args.condition_column, counts_dir)


def run_deseq(counts_matrix: Path, metadata_path: Path, condition_col: str, out_dir: Path):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    counts_df = pd.read_csv(counts_matrix, sep="\t")
    meta_df = pd.read_csv(metadata_path, sep=None, engine="python")
    if "sample" not in meta_df.columns:
        raise ValueError("Metadata must include a 'sample' column.")
    if condition_col not in meta_df.columns:
        raise ValueError(f"Metadata missing condition column: {condition_col}")

    counts_df = counts_df.set_index("gene_id").T
    counts_df.index.name = "sample"
    counts_df = counts_df.astype(int)

    common = counts_df.index.intersection(meta_df["sample"])
    if len(common) < 2:
        raise ValueError("No overlapping samples between counts and metadata.")
    counts_df = counts_df.loc[common]
    meta_df = meta_df.set_index("sample").loc[common]

    dds = DeseqDataSet(counts=counts_df, metadata=meta_df, design_factors=condition_col)
    dds.deseq2()
    stats = DeseqStats(dds, n_cpus=1)
    stats.summary()
    res = stats.results_df.reset_index().rename(columns={"index": "gene_id"})
    out = out_dir / "deseq2_results.tsv"
    res.to_csv(out, sep="\t", index=False)
    print(f"Differential expression results written: {out}")


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
