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

# Citations
If you use MycoGene in your work, please consider this citing repository https://github.com/CDCgov/MycoGene

# CDCgov GitHub Organization Open Source Project

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 

## Access Request, Repo Creation Request

* [CDC GitHub Open Project Request Form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUNk43NzMwODJTRzA4NFpCUk1RRU83RTFNVi4u) _[Requires a CDC Office365 login, if you do not have a CDC Office365 please ask a friend who does to submit the request on your behalf. If you're looking for access to the CDCEnt private organization, please use the [GitHub Enterprise Cloud Access Request form](https://forms.office.com/Pages/ResponsePage.aspx?id=aQjnnNtg_USr6NJ2cHf8j44WSiOI6uNOvdWse4I-C2NUQjVJVDlKS1c0SlhQSUxLNVBaOEZCNUczVS4u).]_

## Related documents

* [Open Practices](open_practices.md)
* [Rules of Behavior](rules_of_behavior.md)
* [Thanks and Acknowledgements](thanks.md)
* [Disclaimer](DISCLAIMER.md)
* [Contribution Notice](CONTRIBUTING.md)
* [Code of Conduct](code-of-conduct.md)

## Overview

Describe the purpose of your project. Add additional sections as necessary to help collaborators and potential collaborators understand and use your project.
  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](DISCLAIMER.md)
and [Code of Conduct](code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).

## SHARE IT Act Metadata
* Organization: NCEZID/DFWED/MDB
* Contact Email: ncezid_shareit@cdc.gov
