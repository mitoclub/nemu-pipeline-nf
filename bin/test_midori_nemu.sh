#!/usr/bin/env bash
# Translate one sequence from a MIDORI BLAST database and run NeMu on it.
#
# The database is any prefix produced by prepare_midori_blastdb.sh. Pass the
# genetic code for that gene and, when the FASTA mixes clades, a lineage
# label so the translation uses the right code. Vertebrate mitochondrial
# genes use --gencode 2 --lineage Mammalia_40674.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: test_midori_nemu.sh --db PREFIX [options]

Options:
  --db PREFIX        BLAST database prefix from prepare_midori_blastdb.sh
  --taxdump DIR      TaxonKit directory (default: $TAXONKIT_DB or ~/.taxonkit)
  --gencode INT      NCBI translation table (default: 2)
  --lineage TEXT     Keep sequences whose FASTA title contains this text
  --min-length INT   Minimum nucleotide length (default: 150)
  --min-seqs INT     Minimum sequences of that species in the database (default: 10)
  --max-seqs INT     Maximum sequences of that species (default: 25)
  --min-hits INT     Species hits required before NeMu starts (default: 4)
  --max-tries INT    Candidates to test before giving up (default: 20)
  --threads INT      CPUs for the NeMu run (default: 4)
  --test-dir DIR     NeMu -output-dir (default: results_midori_test)
  --pipeline FILE    NeMu main.nf (default: this repository's main.nf)
  --seed INT         Shuffle seed for the candidate list
  -h, --help         Show this help

The candidate is a random complete sequence: both gene ends are present, the
length is a multiple of three, the species has a family, and the species has
between --min-seqs and --max-seqs records. A species of that size is enough
for a spectrum and small enough for the alignment to finish.

Example for a vertebrate mitochondrial gene:
  bin/test_midori_nemu.sh \
    --db nuc_db/MIDORI2_UNIQ_NUC_SP_GB274_ND2_NEMU \
    --gencode 2 \
    --lineage Mammalia_40674 \
    --taxdump "$HOME/.taxonkit"
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

DB=""
TAXDUMP="${TAXONKIT_DB:-${HOME}/.taxonkit}"
GENCODE=2
LINEAGE=""
MIN_LENGTH=150
MIN_SEQS=10
MAX_SEQS=25
MIN_HITS=4
MAX_TRIES=20
THREADS=4
TEST_DIR=""
PIPELINE=""
SEED=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --db) DB=$2; shift 2 ;;
        --taxdump) TAXDUMP=$2; shift 2 ;;
        --gencode) GENCODE=$2; shift 2 ;;
        --lineage) LINEAGE=$2; shift 2 ;;
        --min-length) MIN_LENGTH=$2; shift 2 ;;
        --min-seqs) MIN_SEQS=$2; shift 2 ;;
        --max-seqs) MAX_SEQS=$2; shift 2 ;;
        --min-hits) MIN_HITS=$2; shift 2 ;;
        --max-tries) MAX_TRIES=$2; shift 2 ;;
        --threads) THREADS=$2; shift 2 ;;
        --test-dir) TEST_DIR=$2; shift 2 ;;
        --pipeline) PIPELINE=$2; shift 2 ;;
        --seed) SEED=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) die "unknown argument: $1" ;;
    esac
done

[[ -n "$DB" ]] || die "--db is required"
[[ -f "$TAXDUMP/nodes.dmp" && -f "$TAXDUMP/names.dmp" ]] || die "taxdump directory has no nodes.dmp/names.dmp: $TAXDUMP"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
[[ -n "$PIPELINE" ]] || PIPELINE="$REPO_DIR/main.nf"
[[ -n "$TEST_DIR" ]] || TEST_DIR="$REPO_DIR/results_midori_test"

DB=$(readlink -f "$DB")
TAXDUMP=$(readlink -f "$TAXDUMP")
PIPELINE=$(readlink -f "$PIPELINE")
mkdir -p "$TEST_DIR"
TEST_DIR=$(readlink -f "$TEST_DIR")

CLEAN="${DB}.fa"
RECORDS="${DB}.records.tsv"
[[ -s "$CLEAN" && -s "$RECORDS" ]] || die "missing ${CLEAN} or ${RECORDS}; build the database first"
[[ -f "${DB}.nin" || -f "${DB}.ndb" ]] || die "BLAST database not found at $DB"
[[ -f "$PIPELINE" ]] || die "NeMu workflow not found: $PIPELINE"

require_cmd python3
require_cmd taxonkit
require_cmd seqkit
require_cmd tblastn
require_cmd blastdbcmd
require_cmd nextflow

export TAXONKIT_DB="$TAXDUMP"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

echo "Choosing a test sequence from $DB"
cut -f2 "$RECORDS" | sort -u > "$TMP/taxids.txt"
taxonkit reformat2 -I 1 -f '{family}\t{species}' "$TMP/taxids.txt" > "$TMP/ranks.tsv" 2>"$TMP/taxonkit.err"

seed_args=()
if [[ -n "$SEED" ]]; then
    seed_args=(--seed "$SEED")
fi
"$SCRIPT_DIR/midori_test_query.py" candidates \
    --records "$RECORDS" \
    --ranks "$TMP/ranks.tsv" \
    --fasta "$CLEAN" \
    --lineage "$LINEAGE" \
    --min-length "$MIN_LENGTH" \
    --min-seqs "$MIN_SEQS" \
    --max-seqs "$MAX_SEQS" \
    "${seed_args[@]}" \
    --out "$TMP/candidates.tsv"

query=""
seqid=""
taxid=""
family=""
species=""
nhits=""
tries=0
while IFS=$'\t' read -r seqid taxid _seqlen nseq family species; do
    tries=$((tries + 1))
    if [[ "$tries" -gt "$MAX_TRIES" ]]; then
        break
    fi
    echo "  trying $seqid taxid=$taxid species_seqs=$nseq family=$family ($species)"
    "$SCRIPT_DIR/midori_test_query.py" extract \
        --fasta "$CLEAN" --seqid "$seqid" --out "$TMP/query.fa"
    seqkit translate -T "$GENCODE" -f 1 -M -w 0 "$TMP/query.fa" -o "$TMP/query.aa"
    if ! "$SCRIPT_DIR/midori_test_query.py" protein \
        --translation "$TMP/query.aa" --taxid "$taxid" --out "$TMP/query.faa"; then
        continue
    fi

    printf '%s\n' "$taxid" > "$TMP/species.taxid"
    : > "$TMP/species.hits"
    tblastn -query "$TMP/query.faa" -db "$DB" -db_gencode "$GENCODE" \
        -max_target_seqs 50 -evalue 1e-5 \
        -taxidlist "$TMP/species.taxid" -no_taxid_expansion \
        -outfmt '6 saccver pident length qlen sstart send sframe staxids' \
        -out "$TMP/species.hits" >"$TMP/species.out" 2>"$TMP/species.err" || true
    if grep -q "Taxonomy filtering is not supported" "$TMP/species.err"; then
        die "tblastn -taxidlist is not supported on $DB"
    fi
    if ! nhits=$("$SCRIPT_DIR/midori_test_query.py" hits \
        --hits "$TMP/species.hits" --taxid "$taxid" \
        --min-hits "$MIN_HITS" --coords "$TMP/coords.txt"); then
        echo "    skip: $nhits species hits passed the cutoffs"
        continue
    fi
    echo "    species hits passing NeMu cutoffs: $nhits"

    # A taxid that is not in the database must be rejected.
    echo 2 > "$TMP/absent.taxid"
    set +e
    tblastn -query "$TMP/query.faa" -db "$DB" -db_gencode "$GENCODE" \
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

    blastdbcmd -db "$DB" -entry_batch "$TMP/coords.txt" -outfmt '%f' -out "$TMP/fetched.fa"
    [[ -s "$TMP/fetched.fa" ]] || die "blastdbcmd returned no sequence for $(cat "$TMP/coords.txt")"
    query="$TEST_DIR/query.faa"
    cp "$TMP/query.faa" "$query"
    break
done < "$TMP/candidates.tsv"

[[ -n "$query" ]] || die "no candidate produced a usable protein query"
echo "BLAST check passed"
echo "  sequence: $seqid"
echo "  taxid: $taxid ($species)"
echo "  family: $family"
echo "  species hits: $nhits"
echo "  fetched: $(head -n 1 "$TMP/fetched.fa")"
echo "  query: $query"

echo "Running NeMu"
logfile="$TEST_DIR/nextflow.log"
(
    cd "$REPO_DIR"
    export NXF_ANSI_LOG=false
    nextflow run "$PIPELINE" \
        -output-dir "$TEST_DIR" \
        --inputType protein \
        --input "$query" \
        --gencode "$GENCODE" \
        --db "$DB" \
        --taxdump "$TAXDUMP" \
        --threads "$THREADS"
) 2>&1 | tee "$logfile"

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
