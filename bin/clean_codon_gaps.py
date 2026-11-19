#!/usr/bin/env python3
"""Drop gappy codon columns, then drop gappy sequences.

Site filtering keeps or drops a whole codon (three columns) so a coding
alignment stays in frame. A codon is kept when every column has a gap
fraction at or below the site cutoff. A sequence is kept when its gap
fraction in the trimmed alignment is at or below the sequence cutoff.
Gap characters are '-' and '.'. An incomplete trailing codon is dropped.
"""

import sys


def read_fasta(path):
    records = []
    name = None
    seq = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))
                name, seq = line[1:].rstrip("\n"), []
            else:
                seq.append(line.strip())
        if name is not None:
            records.append((name, "".join(seq)))
    return records


def gap_fraction(columns, n_seqs):
    return sum(1 for base in columns if base in "-.") / n_seqs


def clean_codon_gaps(records, site_cut, seq_cut):
    if not records:
        return []

    n_seqs = len(records)
    aln_len = min(len(seq) for _, seq in records)
    keep = []
    for i in range(0, aln_len - (aln_len % 3), 3):
        if all(
            gap_fraction((seq[i + j] for _, seq in records), n_seqs) <= site_cut
            for j in range(3)
        ):
            keep.extend((i, i + 1, i + 2))

    filtered = []
    for name, seq in records:
        trimmed = "".join(seq[i] for i in keep) if keep else ""
        if not trimmed:
            continue
        frac = sum(1 for base in trimmed if base in "-.") / len(trimmed)
        if frac <= seq_cut:
            filtered.append((name, trimmed))
    return filtered


def main():
    try:
        inp, outp, site_cut, seq_cut = sys.argv[1:]
        site_cut, seq_cut = float(site_cut), float(seq_cut)
    except ValueError:
        print(
            "usage: clean_codon_gaps.py in.fasta out.fasta site_cutoff seq_cutoff",
            file=sys.stderr,
        )
        sys.exit(2)

    records = read_fasta(inp)
    filtered = clean_codon_gaps(records, site_cut, seq_cut)
    with open(outp, "w") as handle:
        for name, seq in filtered:
            handle.write(f">{name}\n{seq}\n")


if __name__ == "__main__":
    main()
