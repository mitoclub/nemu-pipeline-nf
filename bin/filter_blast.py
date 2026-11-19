#!/usr/bin/env python3
"""Filter tblastn hits and write blastdbcmd coordinates.

Reads two BLAST outfmt-6 tables (saccver, pident, length, qlen, gapopen,
sstart, send, evalue, bitscore, sframe). Species hits that pass the identity
and coverage cutoffs are kept in file order. The single best outgroup hit by
bitscore is appended. Coordinates are written for `blastdbcmd -entry_batch`.

The log is written to stdout. The pipeline reads the outgroup accession from
this line, so the wording has to stay as it is:

    Selected Outgroup: <accession>, Ident: <pident>%, Score: <bitscore>
"""

import argparse


def parse_blast(filename, min_ident, min_cov):
    hits = []
    try:
        with open(filename) as handle:
            for line in handle:
                parts = line.strip().split("\t")
                if len(parts) < 10:
                    continue
                sacc, pident = parts[0], float(parts[1])
                length, qlen = float(parts[2]), float(parts[3])
                bitscore, sframe = float(parts[8]), int(parts[9])
                sstart, send = parts[5], parts[6]
                coverage = length / qlen
                if pident >= min_ident and coverage >= min_cov:
                    hits.append(
                        {
                            "sacc": sacc,
                            "sstart": sstart,
                            "send": send,
                            "sframe": sframe,
                            "bitscore": bitscore,
                            "pident": pident,
                        }
                    )
    except FileNotFoundError:
        pass
    return hits


def coord_line(hit):
    strand = "minus" if hit["sframe"] < 0 else "plus"
    if strand == "plus":
        start, end = hit["sstart"], hit["send"]
    else:
        start, end = hit["send"], hit["sstart"]
    return f"{hit['sacc']} {start}-{end} {strand}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("species_hits", help="Species tblastn table")
    parser.add_argument("outgroup_hits", help="Outgroup tblastn table")
    parser.add_argument("coords", help="Output path for blastdbcmd coordinates")
    parser.add_argument("--min-ident-species", type=float, default=80.0)
    parser.add_argument("--min-cov-species", type=float, default=0.5)
    parser.add_argument("--min-ident-outgroup", type=float, default=70.0)
    parser.add_argument("--min-cov-outgroup", type=float, default=0.5)
    args = parser.parse_args()

    print("INFO: Filtering BLAST hits")
    hits = parse_blast(args.species_hits, args.min_ident_species, args.min_cov_species)
    outgroup_hits = parse_blast(
        args.outgroup_hits, args.min_ident_outgroup, args.min_cov_outgroup
    )

    if outgroup_hits:
        outgroup_hits.sort(key=lambda hit: hit["bitscore"], reverse=True)
        best_out = outgroup_hits[0]
        print(
            f"Selected Outgroup: {best_out['sacc']}, "
            f"Ident: {best_out['pident']}%, Score: {best_out['bitscore']}"
        )
        hits.append(best_out)
    else:
        print("WARNING: No valid outgroup found.")

    # Count includes the appended outgroup hit, matching the previous log line.
    print(f"Found {len(hits)} valid species hits.")

    with open(args.coords, "w") as handle:
        for hit in hits:
            handle.write(coord_line(hit) + "\n")


if __name__ == "__main__":
    main()
