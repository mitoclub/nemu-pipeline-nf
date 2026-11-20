#!/usr/bin/env python3
"""Print the encoded id of the FASTA record that matches an outgroup id.

Records are paired by order between the original FASTA and the encoded FASTA
(seq_1, seq_2, ...). Matching tries the record id, then whitespace-delimited
tokens, then a substring of the original header. The comparison is
case-sensitive: the id or species name must be passed as it appears in the
header. Prints an empty line when nothing matches. Does not assume the
outgroup is the last record.
"""

import sys


def fasta_headers(path):
    headers = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                headers.append(line[1:].rstrip("\n"))
    return headers


def rec_id(header):
    return header.split()[0] if header else ""


def tokens(header):
    return header.replace("|", " ").replace(";", " ").replace(":", " ").split()


def match_outgroup(orig_headers, enc_ids, outgroup_id):
    token_hit = ""
    substr_hit = ""
    for enc, orig in zip(enc_ids, orig_headers):
        if rec_id(orig) == outgroup_id:
            return enc
        if not token_hit and any(token == outgroup_id for token in tokens(orig)):
            token_hit = enc
        elif not token_hit and not substr_hit and outgroup_id in orig:
            substr_hit = enc
    return token_hit or substr_hit


def main():
    try:
        orig_path, enc_path, outgroup_id = sys.argv[1:]
    except ValueError:
        print(
            "usage: match_outgroup.py original.fasta encoded.fasta OUTGROUP_ID",
            file=sys.stderr,
        )
        sys.exit(2)

    orig_headers = fasta_headers(orig_path)
    enc_ids = [rec_id(header) for header in fasta_headers(enc_path)]
    print(match_outgroup(orig_headers, enc_ids, outgroup_id))


if __name__ == "__main__":
    main()
