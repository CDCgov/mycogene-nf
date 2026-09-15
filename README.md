# 🍄🧬 MycoGene-nf: Gene-Targeted Variant Analysis of Fungal Isolates

![nextflow](https://img.shields.io/badge/nextflow-DSL2-23aa62.svg)
![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)
![platform](https://img.shields.io/badge/platform-Illumina%20%7C%20ONT-blue.svg)

## Introduction

MycoGene-nf is a Nextflow pipeline for identifying nucleotide and amino acid
variants in target genes from fungal whole-genome sequencing data, generated
on either Illumina or Oxford Nanopore Technologies (ONT) platforms.

Given raw reads and one or more reference gene sequences, the pipeline
assembles each sample's genome, locates the target gene(s) by BLAST,
extracts and aligns the matching protein sequence against the reference,
and reports the resulting amino acid mutations — with per-sample QC gating
every step downstream of assembly.

It is built for genes with biologically relevant markers in fungal
pathogens — for example, *CYP51A*, *CYP51B*, *HMG1*, and *hapE* mutations
associated with azole resistance in *Aspergillus fumigatus* — but works for
any gene given a reference protein sequence.

## Summary

- Accepts local FASTQ files, SRA accessions, or both in the same run
- Supports Illumina (paired-end) and ONT (long-read) sequencing data
- Assembles each sample once, regardless of how many genes are queried
- Supports a **single gene** or **multiple genes in one run**, each reported separately
- Per-sample, per-gene QC gate (≥90% reference coverage) before mutations are called
- Optional *Cyp51*-specific promoter tandem-repeat (TR34/TR46/TR53) analysis

## Quick Start

**Requirements:** [Nextflow](https://www.nextflow.io/) (DSL2) and
[Singularity](https://sylabs.io/singularity/).

```bash
nextflow run CDCgov/mycogene-nf \
  --input samplesheet.csv \
  --query_aa cyp51A.fasta \
  --platform illumina \
  --outdir results \
  -profile singularity
```

## Usage

### Reads: choose one or both

| Flag | Description |
|---|---|
| `--input` | CSV samplesheet of local reads (format below) |
| `--add_sra_file` | Plain text file of SRA accessions, one per line — reads are downloaded automatically and combined with `--input` samples if both are given |

At least one of `--input` or `--add_sra_file` is required.

**`--input` samplesheet format** — CSV with a header row:

Illumina (`--platform illumina`):
```csv
sample,fastq_1,fastq_2
SampleA,/path/to/SampleA_R1.fastq.gz,/path/to/SampleA_R2.fastq.gz
SampleB,/path/to/SampleB_R1.fastq.gz,/path/to/SampleB_R2.fastq.gz
```

ONT (`--platform ont`):
```csv
sample,folder
SampleA,/path/to/SampleA_fastq_dir
SampleB,/path/to/SampleB_fastq_dir
```
(`folder` should contain one or more `*.fastq.gz` files for that sample.)

**`--add_sra_file` format** — plain text, one accession per line:
```
SRR25451365
SRR25451359
SRR25451366
```

### Gene query: choose exactly one

| Flag | Description |
|---|---|
| `--query_aa` | Amino acid FASTA for a **single** gene of interest |
| `--multi_query` | CSV for **multiple** genes in one run (format below) |

**`--multi_query` format** — CSV, no header row:
```csv
cyp51A,/path/to/cyp51A.fasta
cyp51B,/path/to/cyp51B.fasta
hmg1,/path/to/hmg1.fasta
hapE,/path/to/hapE.fasta
```
Every sample is BLASTed against every listed gene in the same run. Each
gene's protein FASTA is used as-is (no nucleotide-to-protein translation is
performed — provide amino acid sequences).

### Optional: Cyp51 TR analysis

`--cyp51` runs an additional promoter tandem-repeat analysis
(TR34/TR46/TR53) specific to *CYP51A*. It requires `--query_fa` (the
nucleotide sequence of the gene) and currently **cannot be combined with
`--multi_query`** — run it as a separate single-gene invocation.

## Parameters

| Parameter | Required | Description |
|---|---|---|
| `--input` | one of `--input`/`--add_sra_file` | Samplesheet of local reads (CSV, see above) |
| `--add_sra_file` | one of `--input`/`--add_sra_file` | Text file of SRA accessions to download and include |
| `--platform` | Yes | `illumina` or `ont` |
| `--query_aa` | one of `--query_aa`/`--multi_query` | Amino acid FASTA, single-gene mode |
| `--multi_query` | one of `--query_aa`/`--multi_query` | Gene CSV, multi-gene mode |
| `--outdir` | No (default `results`) | Output directory |
| `--cyp51` | No (default `false`) | Run the Cyp51 TR analysis (single-gene mode only) |
| `--query_fa` | Required if `--cyp51` is set | Nucleotide FASTA for the Cyp51 TR analysis |
| `--help` | No | Print usage and exit |

## Output

```
results/
├── filtered_reads/                    QC-trimmed reads (fastp/fastplong)
├── raw_reads/                         Downloaded SRA reads, if --add_sra_file used
├── Assemblies/<sample>/               Genome assembly (SPAdes/Flye)
├── blast_gene_hits/                   tblastn hits and extracted protein sequences
│   └── <gene>/<sample>/               (gene subfolder only in multi-gene mode)
├── qc_report/qc_report.csv            One row per sample (per gene, in multi-gene mode)
├── plots/                             Protein alignment + SNP visualization (HTML)
│   └── <gene>/                        (gene subfolder only in multi-gene mode)
├── mutation_report/
│   ├── mutations.csv                  Single-gene mode
│   └── <gene>_mutations.csv           One file per gene, multi-gene mode
├── intermediate_outputs/              Intermediate files (e.g. Cyp51 coding/non-coding extracts)
└── Cyp51_analysis/                    Only present if --cyp51 is set
    ├── TR_report/TR_report.tsv        TR34/TR46/TR53/no_TR/unknown_indel call per sample
    └── TR_plot/                       Promoter alignment visualization (HTML)
```

`qc_report.csv` includes a `gene` column (blank in single-gene mode) and a
`qc_status` of `PASS`/`FAIL` per row — only `PASS` rows proceed to
alignment and mutation calling for that sample/gene.

## Run Commands

**Single gene, Illumina, local reads:**
```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --query_aa cyp51A.fasta \
  --platform illumina \
  --outdir results \
  -profile singularity
```

**Single gene, ONT:**
```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --query_aa cyp51A.fasta \
  --platform ont \
  --outdir results \
  -profile singularity
```

**Reads from SRA only:**
```bash
nextflow run main.nf \
  --add_sra_file accessions.txt \
  --query_aa cyp51A.fasta \
  --platform illumina \
  --outdir results \
  -profile singularity
```

**Local reads + SRA reads combined:**
```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --add_sra_file accessions.txt \
  --query_aa cyp51A.fasta \
  --platform illumina \
  --outdir results \
  -profile singularity
```

**Multiple genes in one run:**
```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --multi_query genes.csv \
  --platform illumina \
  --outdir results \
  -profile singularity
```

**With Cyp51 TR analysis (single-gene mode only):**
```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --query_aa cyp51A.fasta \
  --query_fa cyp51A_nt.fasta \
  --cyp51 \
  --platform illumina \
  --outdir results \
  -profile singularity
```

**View all options:**
```bash
nextflow run main.nf --help
```
