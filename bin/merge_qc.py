#!/usr/bin/env python3
"""
merge_qc.py
-----------
Merge per-sample QC CSV files into a single run-level QC report.

Usage:
    merge_qc.py --input sample1_qc.csv sample2_qc.csv ... --output qc_report.csv
"""

import argparse
import csv
import sys

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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",  nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = []
    for csv_file in sorted(args.input):
        with open(csv_file) as fh:
            for row in csv.DictReader(fh):
                rows.append(row)

    if not rows:
        sys.exit("ERROR: no QC rows found in input files")

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    n_fail = sum(1 for r in rows if r.get("qc_status") == "FAIL")
    print(f"QC report written: {args.output} ({len(rows)} samples, {n_fail} failed)", file=sys.stderr)


if __name__ == "__main__":
    main()
