#!/usr/bin/env python3
"""Rewrite a MIDORI unique-nucleotide FASTA for makeblastdb.

Titles look like:

    ACCESSION###root_1;...;Genus_species_TAXID

The accession becomes the BLAST local id. Characters outside the local-id
alphabet are replaced, and ids longer than 50 characters are replaced with
m<n>, because that is the BLAST limit. The original title is kept after the
id. The species taxid is the number at the end of the last lineage field.
"""

import argparse
import re
import sys

TAXID_RE = re.compile(r"_(\d+)$")
UNSAFE_ID_RE = re.compile(r"[^A-Za-z0-9._-]")
MAX_LOCAL_ID = 50


def read_fasta(path):
    header = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def species_taxid(header):
    accession, separator, lineage = header.partition("###")
    if not separator:
        return accession, ""
    last = lineage.split(";")[-1]
    match = TAXID_RE.search(last)
    if not match:
        return accession, ""
    return accession, match.group(1)


def local_id(accession, used, index):
    """Return a unique BLAST local id and a note when the accession was not kept."""
    seqid = UNSAFE_ID_RE.sub("_", accession)
    note = None
    if len(seqid) > MAX_LOCAL_ID:
        note = "short"
        seqid = ""
    elif seqid in used:
        note = "collision"
        seqid = ""
    if seqid:
        return seqid, note
    number = index
    candidate = f"m{number}"
    while candidate in used:
        number += 1
        candidate = f"m{number}"
    return candidate, note


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", help="MIDORI nucleotide FASTA")
    parser.add_argument("--out-fasta", required=True, help="FASTA with BLAST local ids")
    parser.add_argument("--records", required=True, help="seqid, taxid, length, complete")
    args = parser.parse_args()

    used = set()
    written = 0
    missing_taxid = 0
    collisions = 0
    shortened = 0
    complete = 0

    with open(args.out_fasta, "w") as fasta, open(args.records, "w") as records:
        for header, sequence in read_fasta(args.fasta):
            if not header or not sequence:
                continue
            accession, taxid = species_taxid(header)
            if not taxid:
                missing_taxid += 1
                continue
            seqid, note = local_id(accession, used, written + 1)
            if note == "short":
                shortened += 1
            elif note == "collision":
                collisions += 1
            used.add(seqid)
            is_complete = "1" if not any(char in accession for char in "<>,") else "0"
            if is_complete == "1":
                complete += 1
            fasta.write(f">{seqid} {header}\n{sequence}\n")
            records.write(f"{seqid}\t{taxid}\t{len(sequence)}\t{is_complete}\n")
            written += 1
            if written % 50000 == 0:
                print(f"  {written} sequences", file=sys.stderr)

    print(
        f"  {written} sequences, {missing_taxid} without a terminal taxid, "
        f"{collisions} id collisions, {shortened} ids longer than {MAX_LOCAL_ID} characters",
        file=sys.stderr,
    )
    print(f"  complete sequences: {complete}", file=sys.stderr)
    if written == 0:
        sys.exit(1)


if __name__ == "__main__":
    main()
