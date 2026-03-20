# QuickTranscriptome

Minimal RNA-seq workflow with:
- Miniconda installer
- Reference/annotation download from the web
- Alignment with `minimap2`
- Read summarization with `featureCounts`
- Optional downstream differential expression using treatment metadata (`pydeseq2`)

## 1) Install

```bash
cd QuickTranscriptome
chmod +x install.sh
./install.sh
source ./.miniconda/etc/profile.d/conda.sh
conda activate quicktranscriptome
```

## 2) Input FASTQ convention

`--reads-dir` should contain files like:
- Paired-end: `sampleA_R1.fastq.gz`, `sampleA_R2.fastq.gz`
- Single-end: `sampleB.fastq.gz`

## 3) Run pipeline (default species = parapsilosis)

```bash
python quicktranscriptome.py run --reads-dir /path/to/fastqs --out-dir results
```

This will:
1. Download default *Candida parapsilosis* reference FASTA and GFF3
2. Build a minimap2 index
3. Align reads to BAM files
4. Run `featureCounts`
5. Write count matrix to `results/counts/counts_matrix.tsv`

## 4) Optional: provide treatment metadata for DE analysis

Create `metadata.tsv`:

```tsv
sample	treatment
sampleA	control
sampleB	drug
sampleC	control
sampleD	drug
```

Run:

```bash
python quicktranscriptome.py run \
  --reads-dir /path/to/fastqs \
  --out-dir results \
  --metadata metadata.tsv \
  --condition-column treatment
```

DE results are written to:
- `results/counts/deseq2_results.tsv`

## 5) Custom species/reference

If species is not preset, pass URLs directly:

```bash
python quicktranscriptome.py run \
  --reads-dir /path/to/fastqs \
  --species my_species \
  --fasta-url "https://example.org/reference.fa.gz" \
  --gff-url "https://example.org/annotation.gff3.gz"
```
