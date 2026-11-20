#!/usr/bin/env python3
"""Summarize a MIDORI BLAST taxid map against a taxonkit taxdump."""

import argparse
from collections import Counter


def read_table(path):
    with open(path) as handle:
        for line in handle:
            if line.strip():
                yield line.rstrip("\n").split("\t")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", required=True, help="seqid, taxid, length, complete")
    parser.add_argument("--lineage", required=True, help="taxonkit lineage -c -r -n output")
    parser.add_argument("--ranks", required=True, help="taxonkit reformat2 {family}\\t{species}")
    parser.add_argument("--counts", required=True, help="raw taxid class counts")
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--database", required=True)
    parser.add_argument("--taxdump", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rank_of = {}
    for parts in read_table(args.lineage):
        if len(parts) >= 5:
            rank_of[parts[0]] = parts[4]

    family_of = {}
    for parts in read_table(args.ranks):
        taxid = parts[0]
        family_of[taxid] = parts[1] if len(parts) > 1 else ""

    per_taxid = Counter()
    complete = 0
    for parts in read_table(args.records):
        per_taxid[parts[1]] += 1
        if len(parts) > 3 and parts[3] == "1":
            complete += 1

    species_rank = 0
    with_family = 0
    sequences_without_family = 0
    for taxid, count in per_taxid.items():
        if rank_of.get(taxid) == "species":
            species_rank += 1
        if family_of.get(taxid):
            with_family += 1
        else:
            sequences_without_family += count

    lines = [
        "MIDORI BLAST database taxid check",
        f"fasta: {args.fasta}",
        f"database: {args.database}",
        f"taxdump: {args.taxdump}",
        "",
        f"sequences\t{sum(per_taxid.values())}",
        f"unique_taxids\t{len(per_taxid)}",
        f"taxids_ranked_species\t{species_rank}",
        f"taxids_with_family\t{with_family}",
        f"sequences_whose_taxid_has_no_family\t{sequences_without_family}",
        f"complete_sequences\t{complete}",
    ]
    with open(args.counts) as handle:
        lines.extend(line.rstrip("\n") for line in handle if line.strip())

    text = "\n".join(lines) + "\n"
    with open(args.out, "w") as handle:
        handle.write(text)
    print(text, end="")


if __name__ == "__main__":
    main()
