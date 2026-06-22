#!/usr/bin/env python3
"""
sample_qc.py
------------
Parse fastp (Illumina) or fastplong (ONT) JSON, assembly FASTA, and
gene coverage to produce one CSV row of QC metrics per sample.

Primary QC columns appear first; diagnostic columns are pushed to the end.
Also computes an overall qc_status (PASS/FAIL) based on configurable thresholds.

Usage:
    sample_qc.py --sample <id> --platform <illumina|ont> \
                 --json <sample_qc.json> --assembly <scaffolds.fasta> \
                 --gene-coverage <coverage.csv> \
                 --output <sample_qc.csv>
"""

import argparse
import csv
import json
import sys
from pathlib import Path

# ── QC thresholds (fail if any condition is true) ──────────────────────────
MIN_GENE_COVERAGE   = 90.0    # percent
MIN_PERCENT_RETAINED = 80.0   # percent
MIN_N50             = 10000   # bp
MAX_CONTIG_COUNT    = 500     # contigs


def mean_quality_from_curves(data, platform):
    """Average per-base mean quality across the quality_curves section."""
    try:
        if platform == "ont":
            means = data.get("read_quality_curves", {}).get("mean_qual", [])
        else:
            r1 = data.get("read1_quality_curves", {}).get("mean_qual", [])
            r2 = data.get("read2_quality_curves", {}).get("mean_qual", [])
            means = r1 + r2
        if means:
            return round(sum(means) / len(means), 2)
    except Exception:
        pass
    return "NA"


def parse_fastp_json(json_path, platform):
    with open(json_path) as fh:
        data = json.load(fh)

    summary = data.get("summary", {})
    before  = summary.get("before_filtering", {})
    after   = summary.get("after_filtering", {})

    input_reads    = before.get("total_reads", "NA")
    filtered_reads = after.get("total_reads", "NA")
    total_bases    = after.get("total_bases", "NA")

    avg_qscore = mean_quality_from_curves(data, platform)

    gc         = after.get("gc_content", "NA")
    gc_content = round(gc * 100, 2) if gc != "NA" else "NA"

    avg_read_len = after.get("read1_mean_length") or after.get("read_mean_length", "NA")

    if input_reads != "NA" and filtered_reads != "NA" and input_reads > 0:
        pct_retained = round(filtered_reads / input_reads * 100, 2)
    else:
        pct_retained = "NA"

    return {
        "input_reads":      input_reads,
        "filtered_reads":   filtered_reads,
        "percent_retained": pct_retained,
        "avg_qscore":       avg_qscore,
        "gc_content":       gc_content,
        "avg_read_length":  avg_read_len,
        "total_bases":      total_bases,
    }


def parse_assembly_fasta(fasta_path):
    contig_lengths = []
    current_len = 0

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_len > 0:
                    contig_lengths.append(current_len)
                current_len = 0
            else:
                current_len += len(line)
    if current_len > 0:
        contig_lengths.append(current_len)

    if not contig_lengths:
        return {"assembly_size": 0, "contig_count": 0, "N50": 0}

    assembly_size = sum(contig_lengths)
    contig_count  = len(contig_lengths)

    cumulative = 0
    n50 = 0
    for length in sorted(contig_lengths, reverse=True):
        cumulative += length
        if cumulative >= assembly_size / 2:
            n50 = length
            break

    return {
        "assembly_size": assembly_size,
        "contig_count":  contig_count,
        "N50":           n50,
    }


def parse_gene_coverage(coverage_csv):
    with open(coverage_csv) as fh:
        row = next(csv.DictReader(fh))
        return float(row["gene_coverage"])


def determine_qc_status(gene_coverage, percent_retained, n50, contig_count):
    """Return PASS or FAIL with reasons, based on configured thresholds."""
    reasons = []

    if isinstance(gene_coverage, (int, float)) and gene_coverage < MIN_GENE_COVERAGE:
        reasons.append(f"low_gene_coverage(<{MIN_GENE_COVERAGE}%)")
    if isinstance(percent_retained, (int, float)) and percent_retained < MIN_PERCENT_RETAINED:
        reasons.append(f"low_reads_retained(<{MIN_PERCENT_RETAINED}%)")
    if isinstance(n50, (int, float)) and n50 < MIN_N50:
        reasons.append(f"low_N50(<{MIN_N50}bp)")
    if isinstance(contig_count, (int, float)) and contig_count > MAX_CONTIG_COUNT:
        reasons.append(f"high_contig_count(>{MAX_CONTIG_COUNT})")

    if reasons:
        return "FAIL", ";".join(reasons)
    return "PASS", ""


# Primary columns first (at-a-glance), diagnostic columns pushed to the end
FIELDNAMES = [
    "sample",
    "qc_status",
    "qc_fail_reason",
    "gene_coverage",
    "percent_retained",
    "avg_qscore",
    "gc_content",
    "input_reads",
    "filtered_reads",
    # diagnostic / deprioritized columns
    "platform",
    "avg_read_length",
    "total_bases",
    "assembly_size",
    "contig_count",
    "N50",
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample",        required=True)
    parser.add_argument("--platform",      required=True, choices=["illumina", "ont"])
    parser.add_argument("--json",          required=True, help="fastp or fastplong JSON report")
    parser.add_argument("--assembly",      required=True, help="Assembly FASTA")
    parser.add_argument("--gene-coverage", required=True, help="CSV from gene_coverage.py")
    parser.add_argument("--output",        required=True, help="Output CSV path")
    args = parser.parse_args()

    for f in [args.json, args.assembly, args.gene_coverage]:
        if not Path(f).exists():
            sys.exit(f"ERROR: file not found: {f}")

    qc_metrics       = parse_fastp_json(args.json, args.platform)
    assembly_metrics = parse_assembly_fasta(args.assembly)
    gene_coverage    = parse_gene_coverage(args.gene_coverage)

    qc_status, qc_fail_reason = determine_qc_status(
        gene_coverage,
        qc_metrics["percent_retained"],
        assembly_metrics["N50"],
        assembly_metrics["contig_count"],
    )

    row = {
        "sample":         args.sample,
        "platform":       args.platform,
        "qc_status":      qc_status,
        "qc_fail_reason": qc_fail_reason,
        "gene_coverage":  gene_coverage,
        **qc_metrics,
        **assembly_metrics,
    }

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
