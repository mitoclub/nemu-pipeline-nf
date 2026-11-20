#!/usr/bin/env python3
"""Choose and check one nucleotide sequence for a NeMu protein run.

Subcommands:

    candidates   list complete sequences whose species is a usable test
    extract      write one nucleotide record from the cleaned FASTA
    protein      write a NeMu protein query from a seqkit translation
    hits         keep tblastn hits that pass NeMu's species cutoffs
"""

import argparse
import random
import sys
from collections import Counter


def read_tsv(path):
    with open(path) as handle:
        for line in handle:
            if line.strip():
                yield line.rstrip("\n").split("\t")


def command_candidates(args):
    """Write seqid, taxid, length, nseq, family, species."""
    family_of = {}
    species_of = {}
    for parts in read_tsv(args.ranks):
        taxid = parts[0]
        family_of[taxid] = parts[1] if len(parts) > 1 else ""
        species_of[taxid] = parts[2] if len(parts) > 2 else ""

    per_taxid = Counter()
    rows = []
    for parts in read_tsv(args.records):
        seqid, taxid, length = parts[0], parts[1], int(parts[2])
        per_taxid[taxid] += 1
        complete = len(parts) > 3 and parts[3] in {"1", "mammal_complete"}
        rows.append((seqid, taxid, length, complete))

    eligible = {}
    for seqid, taxid, length, complete in rows:
        size = per_taxid[taxid]
        if not complete or length < args.min_length or length % 3:
            continue
        if not family_of.get(taxid):
            continue
        if not (args.min_seqs <= size <= args.max_seqs):
            continue
        eligible[seqid] = (taxid, length, size, family_of[taxid], species_of.get(taxid, ""))

    if args.lineage:
        matched = set()
        prefix = ">"
        with open(args.fasta) as handle:
            for line in handle:
                if not line.startswith(prefix):
                    continue
                seqid, _separator, header = line[1:].rstrip("\n").partition(" ")
                if seqid in eligible and args.lineage in header:
                    matched.add(seqid)
        eligible = {seqid: row for seqid, row in eligible.items() if seqid in matched}

    chosen = list(eligible)
    random.Random(args.seed).shuffle(chosen)
    with open(args.out, "w") as handle:
        for seqid in chosen:
            taxid, length, size, family, species = eligible[seqid]
            handle.write(f"{seqid}\t{taxid}\t{length}\t{size}\t{family}\t{species}\n")
    print(f"  {len(chosen)} candidate sequences", file=sys.stderr)
    if not chosen:
        print("no candidate sequences matched the filters", file=sys.stderr)
        sys.exit(1)


def command_extract(args):
    wanted = f">{args.seqid} "
    writing = False
    with open(args.fasta) as source, open(args.out, "w") as output:
        for line in source:
            if line.startswith(">"):
                if writing:
                    break
                writing = line.startswith(wanted)
            if writing:
                output.write(line)
    if not writing:
        print(f"sequence id not found: {args.seqid}", file=sys.stderr)
        sys.exit(1)


def command_protein(args):
    amino_acids = []
    with open(args.translation) as handle:
        for line in handle:
            if not line.startswith(">"):
                amino_acids.append(line.strip())
    protein = "".join(amino_acids)
    if protein.endswith("*"):
        protein = protein[:-1]
    if len(protein) < 30 or "*" in protein:
        print(
            "skip: translation is shorter than 30 aa or has an internal stop",
            file=sys.stderr,
        )
        sys.exit(1)
    with open(args.out, "w") as handle:
        handle.write(f">query [{args.taxid}]\n{protein}\n")


def command_hits(args):
    """Count species hits and write one blastdbcmd coordinate line.

    The BLAST table columns are:
    saccver, pident, length, qlen, sstart, send, sframe, staxids.
    """
    kept = None
    count = 0
    with open(args.hits) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            sacc, pident, length, qlen = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
            sstart, send, frame, taxid = parts[4], parts[5], int(parts[6]), parts[7]
            if qlen == 0 or pident < args.min_ident or (length / qlen) < args.min_cov:
                continue
            if taxid != args.taxid:
                continue
            count += 1
            if kept is None:
                if frame < 0:
                    start, end, strand = send, sstart, "minus"
                else:
                    start, end, strand = sstart, send, "plus"
                kept = f"{sacc} {start}-{end} {strand}\n"
    print(count)
    if count < args.min_hits or kept is None:
        sys.exit(1)
    with open(args.coords, "w") as handle:
        handle.write(kept)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    candidates = commands.add_parser("candidates", help="list sequences for a test query")
    candidates.add_argument("--records", required=True)
    candidates.add_argument("--ranks", required=True)
    candidates.add_argument("--fasta", required=True, help="cleaned FASTA, used with --lineage")
    candidates.add_argument("--lineage", default="", help="header text that the sequence must contain")
    candidates.add_argument("--min-length", type=int, default=150)
    candidates.add_argument("--min-seqs", type=int, default=10)
    candidates.add_argument("--max-seqs", type=int, default=25)
    candidates.add_argument("--seed", type=int, default=None)
    candidates.add_argument("--out", required=True)
    candidates.set_defaults(func=command_candidates)

    extract = commands.add_parser("extract", help="write one nucleotide record")
    extract.add_argument("--fasta", required=True)
    extract.add_argument("--seqid", required=True)
    extract.add_argument("--out", required=True)
    extract.set_defaults(func=command_extract)

    protein = commands.add_parser("protein", help="write a NeMu protein query")
    protein.add_argument("--translation", required=True)
    protein.add_argument("--taxid", required=True)
    protein.add_argument("--out", required=True)
    protein.set_defaults(func=command_protein)

    hits = commands.add_parser("hits", help="count tblastn hits and write coordinates")
    hits.add_argument("--hits", required=True)
    hits.add_argument("--taxid", required=True)
    hits.add_argument("--min-ident", type=float, default=80.0)
    hits.add_argument("--min-cov", type=float, default=0.5)
    hits.add_argument("--min-hits", type=int, default=4)
    hits.add_argument("--coords", required=True)
    hits.set_defaults(func=command_hits)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
