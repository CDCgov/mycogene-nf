#!/usr/bin/env python3
"""
parse_mutations.py
------------------
Parse a ClustalO protein MSA FASTA and report one row per sample
with all mutations comma-separated.

First sequence in the alignment is treated as the reference.

Output CSV:
    sample    mutation
    SRR001    G98H,K234L
    SRR002    no_mutations
    SRR003    G54W,ins103AA,del200-202

Usage:
    python3 parse_mutations.py --alignment aln_protein_output.fasta \
                               --output mutations.csv
"""

import argparse
import csv
import sys
from pathlib import Path

def parse_fasta(path):
    """Return list of (id, gapped_sequence) from a FASTA file."""
    records = []
    header, seq_parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts).upper()))
    return records

def mutations_for_sample(ref_aln, sample_aln):
    """
    Walk the alignment columns and return a list of mutation strings
    relative to the ungapped reference (1-based positions).

    Examples:
        substitution  -> "G98H"
        insertion     -> "ins98AA"    (AA inserted after ref position 98)
        deletion      -> "del54-56"   (ref positions 54-56 deleted)
    """
    mutations = []
    ref_pos = 0
    i = 0
    n = len(ref_aln)

    while i < n:
        ref_aa = ref_aln[i]
        qry_aa = sample_aln[i]

        # Advance ref position counter only for non-gap ref columns
        if ref_aa != "-":
            ref_pos += 1

        if ref_aa == qry_aa:
            i += 1
            continue

        # -- Insertion: gap in reference ------------------------------------
        if ref_aa == "-":
            ins_seq = ""
            j = i
            while j < n and ref_aln[j] == "-":
                if sample_aln[j] != "-":
                    ins_seq += sample_aln[j]
                j += 1
            if ins_seq:
                mutations.append(f"ins{ref_pos}{ins_seq}")
            i = j
            continue

        # -- Deletion: gap in sample ----------------------------------------
        if qry_aa == "-":
            del_start = ref_pos
            j = i
            while j < n and sample_aln[j] == "-":
                if ref_aln[j] != "-":
                    ref_pos += 1
                j += 1
            del_end = ref_pos
            if del_start == del_end:
                mutations.append(f"del{del_start}")
            else:
                mutations.append(f"del{del_start}-{del_end}")
            i = j
            continue

        # -- Substitution ---------------------------------------------------
        mutations.append(f"{ref_aa}{ref_pos}{qry_aa}")
        i += 1

    return mutations

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--alignment", required=True,
                        help="ClustalO alignment FASTA. First sequence = reference.")
    parser.add_argument("-o", "--output", required=True,
                        help="Output CSV path.")
    args = parser.parse_args()

    if not Path(args.alignment).exists():
        sys.exit(f"ERROR: file not found: {args.alignment}")

    records = parse_fasta(args.alignment)

    if len(records) < 2:
        sys.exit(f"ERROR: need at least 2 sequences (reference + 1 sample), found {len(records)}")

    ref_id, ref_aln = records[0]

    print(f"Reference : {ref_id} ({len(ref_aln.replace('-',''))} aa)", file=sys.stderr)
    print(f"Samples   : {len(records) - 1}", file=sys.stderr)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter=",")
        writer.writerow(["sample", "mutation"])

        for sample_id, sample_aln in records[1:]:
            muts = mutations_for_sample(ref_aln, sample_aln)
            writer.writerow([
                sample_id,
                ",".join(muts) if muts else "no_mutations"
            ])

    print(f"Done. Output: {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()



