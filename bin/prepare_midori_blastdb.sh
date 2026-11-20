#!/usr/bin/env bash
# Rebuild a MIDORI2 unique-nucleotide FASTA into a BLAST database NeMu can
# search with `tblastn -taxidlist -no_taxid_expansion`.
#
# MIDORI titles are:
#   ACCESSION###root_1;...;Genus_species_TAXID
# The prebuilt MIDORI BLAST volume is database version 4 and has no taxid
# index, so -taxidlist is rejected. This script sanitizes the accession,
# takes the terminal species taxid, rewrites merged taxids to the current
# id with taxonkit, and runs:
#   makeblastdb -parse_seqids -taxid_map
#
# With --test, it translates one random complete mammalian ND2 (table 2)
# and runs the NeMu protein workflow against the new database.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: prepare_midori_blastdb.sh --fasta FILE --out PREFIX [options]

Options:
  --fasta FILE       MIDORI nucleotide FASTA
  --out PREFIX       BLAST database prefix to create (must differ from the FASTA prefix)
  --taxdump DIR      TaxonKit data directory with nodes.dmp and names.dmp
                     (default: $TAXONKIT_DB or ~/.taxonkit)
  --test             Run NeMu on one random mammalian ND2 protein
  --test-dir DIR     NeMu -output-dir (default: results_midori_nd2_test)
  --pipeline FILE    NeMu main.nf (default: sibling of this script's parent)
  --gencode INT      Translation table for the test protein (default: 2)
  --threads INT      CPUs for the NeMu test (default: 4)
  -h, --help         Show this help

The test query is a random Mammalia ND2 with both ends complete, length
divisible by 3 and at least 900 bp, from a species that has 10-25 sequences
and a family rank. Vertebrate mitochondrial code 2 matches that clade.
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
RUN_TEST=0
TEST_DIR=""
PIPELINE=""
GENCODE=2
THREADS=4

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fasta) FASTA=$2; shift 2 ;;
        --out) OUT=$2; shift 2 ;;
        --taxdump) TAXDUMP=$2; shift 2 ;;
        --test) RUN_TEST=1; shift ;;
        --test-dir) TEST_DIR=$2; shift 2 ;;
        --pipeline) PIPELINE=$2; shift 2 ;;
        --gencode) GENCODE=$2; shift 2 ;;
        --threads) THREADS=$2; shift 2 ;;
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
REPO_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
[[ -n "$PIPELINE" ]] || PIPELINE="$REPO_DIR/main.nf"
[[ -n "$TEST_DIR" ]] || TEST_DIR="$REPO_DIR/results_midori_nd2_test"
PIPELINE=$(readlink -f "$PIPELINE")
mkdir -p "$(dirname "$OUT")"
mkdir -p "$TEST_DIR"
TEST_DIR=$(readlink -f "$TEST_DIR")

require_cmd seqkit
require_cmd taxonkit
require_cmd makeblastdb
require_cmd blastdbcmd
require_cmd tblastn
require_cmd awk
require_cmd shuf
if [[ "$RUN_TEST" -eq 1 ]]; then
    require_cmd nextflow
    [[ -f "$PIPELINE" ]] || die "NeMu workflow not found: $PIPELINE"
fi

export TAXONKIT_DB="$TAXDUMP"

CLEAN="${OUT}.fa"
MAP="${OUT}.taxid_map"
RECORDS="${OUT}.records.tsv"
REPORT="${OUT}.taxid_check.txt"
QUERY="${OUT}.test_query.faa"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

echo "Parsing titles and species taxids"
# Column 1 of seqkit fx2tab is the full FASTA header. The accession is the
# token before ###; the species taxid is the trailing _<digits> of the lineage.
seqkit fx2tab "$FASTA" \
| awk -F '\t' -v fa="$CLEAN" -v raw="$TMP/raw.map" -v rec="$TMP/records.raw.tsv" '
BEGIN { notax = 0; dups = 0; nshort = 0; n = 0 }
{
    header = $1
    seq = $2
    if (header == "" || seq == "") next
    acc = header
    lineage = ""
    hash = index(header, "###")
    if (hash > 0) {
        acc = substr(header, 1, hash - 1)
        lineage = substr(header, hash + 3)
    }
    nrank = split(lineage, ranks, ";")
    last = ranks[nrank]
    tax = ""
    if (match(last, /_[0-9]+$/)) {
        tax = substr(last, RSTART + 1, RLENGTH - 1)
    }
    if (tax == "") {
        notax++
        next
    }
    id = acc
    gsub(/[^A-Za-z0-9._-]/, "_", id)
    # BLAST local ids are limited to 50 characters.
    if (length(id) > 50) {
        nshort++
        id = ""
    }
    if (id != "" && (id in used)) {
        dups++
        id = ""
    }
    if (id == "") {
        k = n + 1
        do {
            id = sprintf("m%d", k)
            k++
        } while (id in used)
    }
    used[id] = 1
    tag = "-"
    if (index(header, "Mammalia_40674") > 0 && acc !~ /[<>,]/) {
        tag = "mammal_complete"
    }
    print ">" id " " header > fa
    print seq > fa
    print id, tax > raw
    print id "\t" tax "\t" length(seq) "\t" tag > rec
    n++
    if (n % 50000 == 0) printf "  %d sequences\n", n > "/dev/stderr"
}
END {
    printf "  %d sequences, %d without a terminal taxid, %d id collisions, %d ids longer than 50 characters\n", n, notax, dups, nshort > "/dev/stderr"
    if (n == 0) exit 1
}
'

echo "Checking taxids with taxonkit"
awk '{print $2}' "$TMP/raw.map" | sort -u > "$TMP/raw.taxids"
# -L drops the lineage string. -r is required with -L, and -c replaces
# merged taxids with the current id (or -1 / 0 when the node is gone).
taxonkit lineage -c -r -L "$TMP/raw.taxids" > "$TMP/status.tsv" 2>"$TMP/taxonkit_status.err"

# status column: -1 missing, 0 deleted, >0 current taxid (itself, or the new id of a merged node)
awk 'NR == FNR {
    status[$1] = $2
    next
}
{
    id = $1
    raw = $2
    final = raw
    cls = "missing"
    if (raw in status && status[raw] + 0 > 0) {
        final = status[raw]
        cls = (final == raw ? "current" : "merged")
    } else if (raw in status && status[raw] + 0 == 0) {
        cls = "deleted"
    }
    print id, final
    print raw "\t" final "\t" cls > "'"$TMP/taxid_class.tsv"'"
}' "$TMP/status.tsv" "$TMP/raw.map" > "$MAP"

awk -F '\t' 'NR == FNR {
    final[$1] = $2
    next
}
{
    id = $1
    raw = $2
    if (raw in final) $2 = final[raw]
    print $1 "\t" $2 "\t" $3 "\t" $4
}' "$TMP/taxid_class.tsv" "$TMP/records.raw.tsv" > "$RECORDS"

awk '{print $2}' "$MAP" | sort -u > "$TMP/final.taxids"
taxonkit lineage -c -r -n "$TMP/final.taxids" > "$TMP/lineage.tsv" 2>>"$TMP/taxonkit_status.err"
taxonkit reformat2 -I 1 -f '{family}\t{species}' "$TMP/final.taxids" > "$TMP/ranks.tsv" 2>>"$TMP/taxonkit_status.err"

awk -F '\t' -v lineage="$TMP/lineage.tsv" -v ranks="$TMP/ranks.tsv" '
FILENAME == lineage {
    rank[$1] = $5
    name[$1] = $4
    next
}
FILENAME == ranks {
    family[$1] = $2
    species[$1] = $3
    next
}
{
    tax = $2
    nseq[tax]++
    tag = $4
    ntag[tag]++
}
END {
    seqs = 0
    for (t in nseq) {
        seqs += nseq[t]
        nt++
        if (rank[t] == "species") n_species_rank++
        if (family[t] != "") n_with_family++
        else seq_no_family += nseq[t]
    }
    print "sequences\t" seqs
    print "unique_taxids\t" nt
    print "taxids_ranked_species\t" n_species_rank + 0
    print "taxids_with_family\t" n_with_family + 0
    print "sequences_whose_taxid_has_no_family\t" seq_no_family + 0
    print "mammal_complete_sequences\t" ntag["mammal_complete"] + 0
}
' "$TMP/lineage.tsv" "$TMP/ranks.tsv" "$RECORDS" > "$TMP/summary.tsv"

awk -F '\t' '
$3 == "current" { current[$1] = 1 }
$3 == "merged" { merged[$1] = 1 }
$3 == "deleted" { deleted[$1] = 1 }
$3 == "missing" { missing[$1] = 1 }
END {
    print "raw_taxids_current\t" length(current)
    print "raw_taxids_merged\t" length(merged)
    print "raw_taxids_deleted\t" length(deleted)
    print "raw_taxids_missing\t" length(missing)
}
' "$TMP/taxid_class.tsv" > "$TMP/raw_summary.tsv"

{
    echo "MIDORI BLAST database taxid check"
    echo "fasta: $FASTA"
    echo "database: $OUT"
    echo "taxdump: $TAXDUMP"
    echo
    cat "$TMP/summary.tsv" "$TMP/raw_summary.tsv"
} > "$REPORT"
cat "$REPORT"

echo "Building BLAST database"
rm -f "${OUT}".nhr "${OUT}".nin "${OUT}".nsq "${OUT}".njs \
    "${OUT}".ndb "${OUT}".nog "${OUT}".nos "${OUT}".not \
    "${OUT}".ntf "${OUT}".nto "${OUT}".ndb
makeblastdb \
    -in "$CLEAN" \
    -dbtype nucl \
    -parse_seqids \
    -taxid_map "$MAP" \
    -out "$OUT" \
    -title "MIDORI2 UNIQ NUC SP GB274 ND2, species taxids for NeMu" \
    > "$TMP/makeblastdb.log"
cat "$TMP/makeblastdb.log"
grep -q '"number-of-taxids"' "${OUT}.njs" || die "makeblastdb did not store taxids (${OUT}.njs)"
echo "number-of-taxids: $(sed -n 's/.*"number-of-taxids": \([0-9]*\).*/\1/p' "${OUT}.njs")"

echo "Sampling a mammalian ND2 for the BLAST check"
# Species size is the full taxid count. NeMu will retrieve every sequence
# of that taxid, so a species with tens of thousands of records is a poor test.
awk -F '\t' '{c[$2]++} END {for (t in c) print t "\t" c[t]}' "$RECORDS" > "$TMP/counts.tsv"
awk -F '\t' 'NR == FNR {fam[$1] = $2; next} {print $0 "\t" fam[$2]}' "$TMP/ranks.tsv" "$RECORDS" \
| awk -F '\t' 'NR == FNR {n[$1] = $2; next}
$4 == "mammal_complete" && $3 >= 900 && $3 % 3 == 0 && $5 != "" && n[$2] >= 10 && n[$2] <= 25 {
    print $1 "\t" $2 "\t" $3 "\t" n[$2] "\t" $5
}' "$TMP/counts.tsv" - > "$TMP/candidates.tsv"

n_cand=$(wc -l < "$TMP/candidates.tsv")
[[ "$n_cand" -gt 0 ]] || die "no mammalian ND2 candidate with 10-25 sequences and a family rank"
echo "  $n_cand candidate sequences"
shuf "$TMP/candidates.tsv" > "$TMP/candidates.shuffled.tsv"

picked=""
while IFS=$'\t' read -r seqid taxid seqlen nseq family; do
    echo "  trying $seqid taxid=$taxid length=$seqlen species_seqs=$nseq family=$family"
    pat=$(printf '%s\n' "$seqid" | sed -e 's/[][\.^$*+?{}|()]/\\&/g')
    seqkit grep -r -p "^${pat}$" "$CLEAN" > "$TMP/query.fa" || true
    [[ -s "$TMP/query.fa" ]] || continue
    seqkit translate -T "$GENCODE" -f 1 -M -w 0 "$TMP/query.fa" | awk 'NR > 1' > "$TMP/query.aa"
    prot=$(tr -d '\n' < "$TMP/query.aa")
    prot=${prot%\*}
    if [[ ${#prot} -lt 30 || "$prot" == *'*'* ]]; then
        echo "    skip: translation is shorter than 30 aa or has an internal stop"
        continue
    fi
    printf '>ND2 [%s]\n%s\n' "$taxid" "$prot" > "$TMP/query.faa"
    printf '%s\n' "$taxid" > "$TMP/species.taxid"
    : > "$TMP/species.hits"
    tblastn -query "$TMP/query.faa" -db "$OUT" -db_gencode "$GENCODE" \
        -max_target_seqs 50 -evalue 1e-5 \
        -taxidlist "$TMP/species.taxid" -no_taxid_expansion \
        -outfmt '6 saccver pident length qlen sstart send sframe staxids' \
        -out "$TMP/species.hits" >"$TMP/species.out" 2>"$TMP/species.err" || true
    if grep -q "Taxonomy filtering is not supported" "$TMP/species.err"; then
        die "tblastn -taxidlist is not supported on $OUT"
    fi
    awk -F '\t' -v tax="$taxid" '
        $2 >= 80 && ($3 / $4) >= 0.5 && $8 == tax { n++ }
        END { print n + 0 }
    ' "$TMP/species.hits" > "$TMP/nhits"
    nhits=$(cat "$TMP/nhits")
    echo "    species hits passing NeMu cutoffs: $nhits"
    if [[ "$nhits" -lt 4 ]]; then
        echo "    skip: fewer than 4 species hits"
        continue
    fi
    # A taxid absent from the database must be rejected, which is what
    # -taxidlist does when the index is present.
    echo 2 > "$TMP/absent.taxid"
    set +e
    tblastn -query "$TMP/query.faa" -db "$OUT" -db_gencode "$GENCODE" \
        -max_target_seqs 5 -evalue 1e-5 \
        -taxidlist "$TMP/absent.taxid" -no_taxid_expansion \
        -outfmt '6 saccver staxids' \
        -out "$TMP/absent.hits" >"$TMP/absent.out" 2>"$TMP/absent.err"
    absent_status=$?
    set -e
    if [[ "$absent_status" -eq 0 ]]; then
        die "tblastn -taxidlist accepted taxid 2, which is not in this database"
    fi
    grep -q "Taxonomy ID" "$TMP/absent.err" || die "unexpected taxidlist error: $(cat "$TMP/absent.err")"

    awk -F '\t' -v tax="$taxid" '
        $2 >= 80 && ($3 / $4) >= 0.5 && $8 == tax {
            strand = ($7 < 0 ? "minus" : "plus")
            if ($7 < 0) { start = $6; end = $5 } else { start = $5; end = $6 }
            print $1, start "-" end, strand
            exit
        }
    ' "$TMP/species.hits" > "$TMP/coords.txt"
    blastdbcmd -db "$OUT" -entry_batch "$TMP/coords.txt" -outfmt '%f' -out "$TMP/fetched.fa"
    [[ -s "$TMP/fetched.fa" ]] || die "blastdbcmd returned no sequence for $(cat "$TMP/coords.txt")"
    picked=$seqid
    cp "$TMP/query.faa" "$QUERY"
    break
done < "$TMP/candidates.shuffled.tsv"

[[ -n "$picked" ]] || die "no candidate produced a self-hit under NeMu identity and coverage cutoffs"
name=$(awk -F '\t' -v t="$taxid" '$1 == t { print $4; exit }' "$TMP/lineage.tsv")
echo "BLAST check passed"
echo "  sequence: $picked"
echo "  taxid: $taxid ($name)"
echo "  family: $family"
echo "  species hits: $nhits"
echo "  fetched: $(head -n 1 "$TMP/fetched.fa")"
echo "  query: $QUERY"

if [[ "$RUN_TEST" -eq 0 ]]; then
    echo "Database is ready: $OUT"
    echo "Re-run with --test to execute NeMu on $QUERY"
    exit 0
fi

echo "Running NeMu"
echo "  pipeline: $PIPELINE"
echo "  output: $TEST_DIR"
logfile="${OUT}.nemu.log"
(
    cd "$REPO_DIR"
    export NXF_ANSI_LOG=false
    nextflow run "$PIPELINE" \
        -output-dir "$TEST_DIR" \
        --inputType protein \
        --input "$QUERY" \
        --gencode "$GENCODE" \
        --db "$OUT" \
        --taxdump "$TAXDUMP" \
        --threads "$THREADS"
) 2>&1 | tee "$logfile"

# The outgroup line is written by the BLAST task, not by the Nextflow console.
if [[ ! -s "$TEST_DIR/spectra_total.tsv" ]]; then
    die "NeMu did not write $TEST_DIR/spectra_total.tsv. See $logfile"
fi
outgroup_log=$(find "$REPO_DIR/work" -name blast_filtering_log.txt -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -1 | cut -d' ' -f2-)
if [[ -n "$outgroup_log" ]] && grep -q "Selected Outgroup:" "$outgroup_log"; then
    echo "NeMu selected an outgroup:"
    grep "Selected Outgroup:" "$outgroup_log"
else
    echo "WARNING: spectra were written, but no outgroup line was found in the task log."
fi
echo "NeMu test output: $TEST_DIR"
