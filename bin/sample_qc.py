#!/usr/bin/env python3
"""
sample_qc.py
------------
Parse fastp (Illumina) or fastplong (ONT) JSON, assembly FASTA, and
gene coverage to produce one CSV row of QC metrics per sample.

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

# ── QC thresholds ──────────────────────────────────────────────────────────
MIN_GENE_COVERAGE = 90.0   # percent
MIN_AVG_QSCORE    = 20.0   # mean quality score

FIELDNAMES = [
    "sample",
    "gene_coverage",
    "avg_qscore",
    "percent_retained",
    "gc_content",
    "avg_read_length",
    "assembly_size",
    "qc_status",
    "qc_fail_reason",
]


def mean_quality_from_curves(data, platform):
    try:
        if platform == "ont":
            means = data.get("read_quality_curves", {}).get("mean_qual", [])
        else:
            r1 = data.get("read1_before_filtering", {}).get("quality_curves", {}).get("mean", [])
            r2 = data.get("read2_befire_filtering", {}).get("quality_curves", {}).get("mean", [])
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

    gc         = after.get("gc_content", "NA")
    gc_content = round(gc * 100, 2) if gc != "NA" else "NA"

    avg_read_len = after.get("read1_mean_length") or after.get("read_mean_length", "NA")
    avg_qscore   = mean_quality_from_curves(data, platform)

    if input_reads != "NA" and filtered_reads != "NA" and input_reads > 0:
        pct_retained = round(filtered_reads / input_reads * 100, 2)
    else:
        pct_retained = "NA"

    return {
        "percent_retained": pct_retained,
        "avg_qscore":       avg_qscore,
        "gc_content":       gc_content,
        "avg_read_length":  avg_read_len,
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
        return {"assembly_size": 0}

    return {"assembly_size": sum(contig_lengths)}


def parse_gene_coverage(coverage_csv):
    with open(coverage_csv) as fh:
        row = next(csv.DictReader(fh))
        return float(row["gene_coverage"])


def determine_qc_status(gene_coverage, avg_qscore):
    reasons = []

    if isinstance(gene_coverage, (int, float)) and gene_coverage < MIN_GENE_COVERAGE:
        reasons.append(f"low_gene_coverage(<{MIN_GENE_COVERAGE}%)")
    if isinstance(avg_qscore, (int, float)) and avg_qscore < MIN_AVG_QSCORE:
        reasons.append(f"low_avg_qscore(<{MIN_AVG_QSCORE})")

    if reasons:
        return "FAIL", ";".join(reasons)
    return "PASS", ""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample",        required=True)
    parser.add_argument("--platform",      required=True, choices=["illumina", "ont"])
    parser.add_argument("--json",          required=True)
    parser.add_argument("--assembly",      required=True)
    parser.add_argument("--gene-coverage", required=True)
    parser.add_argument("--output",        required=True)
    args = parser.parse_args()

    for f in [args.json, args.assembly, args.gene_coverage]:
        if not Path(f).exists():
            sys.exit(f"ERROR: file not found: {f}")

    qc_metrics       = parse_fastp_json(args.json, args.platform)
    assembly_metrics = parse_assembly_fasta(args.assembly)
    gene_coverage    = parse_gene_coverage(args.gene_coverage)

    qc_status, qc_fail_reason = determine_qc_status(
        gene_coverage,
        qc_metrics["avg_qscore"],
    )

    row = {
        "sample":           args.sample,
        "gene_coverage":    gene_coverage,
        "avg_qscore":       qc_metrics["avg_qscore"],
        "percent_retained": qc_metrics["percent_retained"],
        "gc_content":       qc_metrics["gc_content"],
        "avg_read_length":  qc_metrics["avg_read_length"],
        "assembly_size":    assembly_metrics["assembly_size"],
        "qc_status":        qc_status,
        "qc_fail_reason":   qc_fail_reason,
    }

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
