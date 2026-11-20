#!/usr/bin/env python3
"""Replace MIDORI taxids with the current taxonkit id.

`taxonkit lineage -c -r -L` writes taxid, status, rank. A status above zero
is the current taxid: the same id, or the new id of a merged node. Zero means
the node was deleted and -1 means it is absent. Those taxids are left as they
appear in the FASTA, so the sequence stays in the BLAST database.
"""

import argparse
from collections import Counter


def read_status(path):
    status = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            taxid, code, *_rest = line.rstrip("\n").split("\t")
            status[taxid] = code
    return status


def current_taxid(raw, status):
    code = status.get(raw)
    if code is None:
        return raw, "missing"
    if int(code) > 0:
        kind = "current" if code == raw else "merged"
        return code, kind
    if int(code) == 0:
        return raw, "deleted"
    return raw, "missing"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", required=True, help="seqid, raw taxid, length, complete")
    parser.add_argument("--status", required=True, help="taxonkit lineage -c -r -L output")
    parser.add_argument("--out-records", required=True)
    parser.add_argument("--out-map", required=True, help="seqid and taxid for makeblastdb -taxid_map")
    parser.add_argument("--out-counts", required=True, help="unique raw-taxid class counts")
    args = parser.parse_args()

    status = read_status(args.status)
    classes = Counter()
    seen = {}

    with open(args.records) as source, open(args.out_records, "w") as records, open(args.out_map, "w") as mapping:
        for line in source:
            parts = line.rstrip("\n").split("\t")
            seqid, raw = parts[0], parts[1]
            final, kind = current_taxid(raw, status)
            if raw not in seen:
                seen[raw] = kind
                classes[kind] += 1
            parts[1] = final
            records.write("\t".join(parts) + "\n")
            mapping.write(f"{seqid} {final}\n")

    labels = (
        ("raw_taxids_current", "current"),
        ("raw_taxids_merged", "merged"),
        ("raw_taxids_deleted", "deleted"),
        ("raw_taxids_missing", "missing"),
    )
    with open(args.out_counts, "w") as handle:
        for label, kind in labels:
            handle.write(f"{label}\t{classes[kind]}\n")


if __name__ == "__main__":
    main()
