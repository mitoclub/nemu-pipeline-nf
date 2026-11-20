#!/usr/bin/env bash
# Build a BLAST database NeMu can search with tblastn -taxidlist.
#
# MIDORI unique-nucleotide titles carry a species taxid, but the BLAST
# volumes shipped with MIDORI are database version 4 and have no taxid index.
# This script writes BLASTDB version 5 with -parse_seqids and -taxid_map.
# The FASTA can be any MIDORI gene, not only ND2.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: prepare_midori_blastdb.sh --fasta FILE --out PREFIX [options]

Options:
  --fasta FILE     MIDORI nucleotide FASTA
  --out PREFIX     BLAST database prefix to create
  --taxdump DIR    TaxonKit directory with nodes.dmp and names.dmp
                   (default: $TAXONKIT_DB or ~/.taxonkit)
  -h, --help       Show this help

Example:
  bin/prepare_midori_blastdb.sh \
    --fasta nuc_db/MIDORI2_UNIQ_NUC_SP_GB274_ND2_BLAST.fasta \
    --out nuc_db/MIDORI2_UNIQ_NUC_SP_GB274_ND2_NEMU \
    --taxdump "$HOME/.taxonkit"

The same command works for the other MIDORI gene FASTA files. Run
bin/test_midori_nemu.sh against the new prefix to check it.
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

FASTA=""
OUT=""
TAXDUMP="${TAXONKIT_DB:-${HOME}/.taxonkit}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta) FASTA=$2; shift 2 ;;
        --out) OUT=$2; shift 2 ;;
        --taxdump) TAXDUMP=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown argument: $1" ;;
    esac
done

[[ -n "$FASTA" && -s "$FASTA" ]] || die "--fasta must be an existing FASTA file"
[[ -n "$OUT" ]] || die "--out is required"
[[ -f "$TAXDUMP/nodes.dmp" && -f "$TAXDUMP/names.dmp" ]] || die "taxdump directory has no nodes.dmp/names.dmp: $TAXDUMP"

FASTA=$(readlink -f "$FASTA")
OUT=$(readlink -f "$OUT")
TAXDUMP=$(readlink -f "$TAXDUMP")
src_prefix=${FASTA%.fasta}
src_prefix=${src_prefix%.fa}
src_prefix=${src_prefix%.fna}
[[ "$OUT" != "$src_prefix" ]] || die "--out must not reuse the source FASTA prefix ($src_prefix)"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$(dirname "$OUT")"

require_cmd python3
require_cmd taxonkit
require_cmd makeblastdb
require_cmd cut
require_cmd sort

export TAXONKIT_DB="$TAXDUMP"

CLEAN="${OUT}.fa"
MAP="${OUT}.taxid_map"
RECORDS="${OUT}.records.tsv"
REPORT="${OUT}.taxid_check.txt"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

echo "Parsing titles and species taxids"
"$SCRIPT_DIR/midori_prepare_fasta.py" "$FASTA" \
    --out-fasta "$CLEAN" \
    --records "$TMP/records.raw.tsv"

echo "Checking taxids with taxonkit"
cut -f2 "$TMP/records.raw.tsv" | sort -u > "$TMP/raw.taxids"
# -c rewrites merged taxids. -L omits the lineage string and needs -r.
taxonkit lineage -c -r -L "$TMP/raw.taxids" > "$TMP/status.tsv" 2>"$TMP/taxonkit.err"
"$SCRIPT_DIR/midori_apply_taxids.py" \
    --records "$TMP/records.raw.tsv" \
    --status "$TMP/status.tsv" \
    --out-records "$RECORDS" \
    --out-map "$MAP" \
    --out-counts "$TMP/class_counts.tsv"

cut -d' ' -f2 "$MAP" | sort -u > "$TMP/final.taxids"
taxonkit lineage -c -r -n "$TMP/final.taxids" > "$TMP/lineage.tsv" 2>>"$TMP/taxonkit.err"
taxonkit reformat2 -I 1 -f '{family}\t{species}' "$TMP/final.taxids" > "$TMP/ranks.tsv" 2>>"$TMP/taxonkit.err"
"$SCRIPT_DIR/midori_taxid_report.py" \
    --records "$RECORDS" \
    --lineage "$TMP/lineage.tsv" \
    --ranks "$TMP/ranks.tsv" \
    --counts "$TMP/class_counts.tsv" \
    --fasta "$FASTA" \
    --database "$OUT" \
    --taxdump "$TAXDUMP" \
    --out "$REPORT"

echo "Building BLAST database"
rm -f "${OUT}".nhr "${OUT}".nin "${OUT}".nsq "${OUT}".njs \
    "${OUT}".ndb "${OUT}".nog "${OUT}".nos "${OUT}".not \
    "${OUT}".ntf "${OUT}".nto
makeblastdb \
    -in "$CLEAN" \
    -dbtype nucl \
    -parse_seqids \
    -taxid_map "$MAP" \
    -out "$OUT" \
    -title "$(basename "$OUT")" \
    > "$TMP/makeblastdb.log"
cat "$TMP/makeblastdb.log"
grep -q '"number-of-taxids"' "${OUT}.njs" || die "makeblastdb did not store taxids (${OUT}.njs)"
echo "number-of-taxids: $(sed -n 's/.*"number-of-taxids": \([0-9]*\).*/\1/p' "${OUT}.njs")"
echo "Database is ready: $OUT"
echo "Taxid report: $REPORT"
