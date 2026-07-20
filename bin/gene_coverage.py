#!/usr/bin/env python3
"""
gene_coverage.py
----------------
Calculate identity-weighted gene coverage from tblastn output.

Coverage = (sum of matched positions across all HSPs) / reference_length * 100
where matched positions per HSP = hsp_length * (pident / 100)

Usage:
    gene_coverage.py --sample <id> --tblastn-out <tblastn.tsv> \
                     --ref-length <int> --output <coverage.csv>
"""

import argparse
import csv
import sys
from pathlib import Path


def calculate_coverage(tblastn_tsv, ref_length):
    """
    tblastn_tsv columns expected: qstart qend pident length
    Returns identity-weighted coverage as a percentage of ref_length.
    """
    matched_positions = 0.0

    with open(tblastn_tsv) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            pident   = float(fields[2])
            hsp_len  = float(fields[3])
            matched_positions += hsp_len * (pident / 100.0)

    if ref_length <= 0:
        return 0.0

    coverage = (matched_positions / ref_length) * 100.0
    return round(min(coverage, 100.0), 2)  # cap at 100 in case of overlapping HSPs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample",      required=True)
    parser.add_argument("--tblastn-out", required=True, help="tblastn output: qstart qend pident length")
    parser.add_argument("--ref-length",  required=True, type=int, help="Reference amino acid length")
    parser.add_argument("--output",      required=True)
    args = parser.parse_args()

    if not Path(args.tblastn_out).exists():
        sys.exit(f"ERROR: file not found: {args.tblastn_out}")

    coverage = calculate_coverage(args.tblastn_out, args.ref_length)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "gene_coverage"])
        writer.writerow([args.sample, coverage])


if __name__ == "__main__":
    main()
