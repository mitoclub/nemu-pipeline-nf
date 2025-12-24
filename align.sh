#!/bin/bash

# =================CONFIGURATION=================
# INPUT/OUTPUT
INPUT_FILE="$1"          # First argument: Input FASTA
GENCODE="${2:-1}"        # Second argument: Genetic Code (default: 1)
OUTPUT_PREFIX="final_output"
THREADS=4

# PATHS (Adjust these to your system)
MACSE_JAR="/opt/macse_v2.07.jar"
# Ensure mafft, goalign, and seqkit are in your $PATH

# THRESHOLDS
LARGE_DATA_CUTOFF=500    # Switch to Big Data workflow if seqs > this
MIN_SEQ_LEN=100          # Pre-filter: remove sequences shorter than 100bp
MAX_GAP_SEQ=0.50         # Post-filter: Remove SEQS with >50% gaps
MAX_GAP_SITE=0.50        # Post-filter: Remove SITES (columns) with >50% gaps

# ===============================================

if [[ -z "$INPUT_FILE" ]]; then
    echo "Usage: ./align_and_filter.sh <input.fasta> [gencode]"
    exit 1
fi

echo "--- STARTING PIPELINE ---"
echo "Input: $INPUT_FILE | Gencode: $GENCODE"

# 0. PRE-FILTERING (Sanity Check)
# Remove very short fragments before we waste time aligning them
echo "[0/4] Pre-filtering short sequences (<${MIN_SEQ_LEN}bp)..."
seqkit seq -m $MIN_SEQ_LEN -g "$INPUT_FILE" > input_clean.fasta
SEQ_COUNT=$(grep -c "^>" input_clean.fasta)
echo "      Sequences remaining: $SEQ_COUNT"

# 1. ALIGNMENT LOGIC
if [ "$SEQ_COUNT" -gt "$LARGE_DATA_CUTOFF" ]; then
    # ================= STRATEGY A: BIG DATA WORKFLOW =================
    echo "[1/4] Strategy: BIG DATA (> $LARGE_DATA_CUTOFF sequences)"
    
    # A1. Trim non-homologous fragments
    java -jar "$MACSE_JAR" -prog trimNonHomologousFragments \
        -seq input_clean.fasta -gc_def "$GENCODE" \
        -out_NT 1_trimmed.fasta > /dev/null 2>&1

    # A2. Repair Frameshifts (Fast mode)
    echo "      Running MACSE Repair (Frameshift detection)..."
    java -jar "$MACSE_JAR" -prog alignSequences \
        -seq 1_trimmed.fasta -gc_def "$GENCODE" \
        -out_NT 2_repaired.fasta \
        -max_refine_iter 0 -local_realign_init 0 > /dev/null 2>&1

    # A3. Prepare for Translation (Remove alignment gaps '-', keep FS fixes '!')
    sed 's/-//g' 2_repaired.fasta > 3_ungapped.fasta

    # A4. Translate
    java -jar "$MACSE_JAR" -prog translateNT2AA \
        -seq 3_ungapped.fasta -gc_def "$GENCODE" \
        -out_AA 4_protein.faa > /dev/null 2>&1

    # A5. Align Protein (MAFFT)
    echo "      Running MAFFT Alignment..."
    mafft --thread "$THREADS" --auto --quiet 4_protein.faa > 5_aligned_protein.faa

    # A6. Back-Translate
    java -jar "$MACSE_JAR" -prog reportGapsAA2NT \
        -align_AA 5_aligned_protein.faa \
        -seq 3_ungapped.fasta -gc_def "$GENCODE" \
        -out_NT raw_alignment.fasta > /dev/null 2>&1

    # Cleanup intermediate files
    rm 1_trimmed.fasta 2_repaired.fasta 3_ungapped.fasta 4_protein.faa 5_aligned_protein.faa

else
    # ================= STRATEGY B: PURE MACSE =================
    echo "[1/4] Strategy: PURE MACSE (< $LARGE_DATA_CUTOFF sequences)"
    
    echo "      Running Full MACSE Alignment..."
    java -jar "$MACSE_JAR" -prog alignSequences \
        -seq input_clean.fasta -gc_def "$GENCODE" \
        -out_NT raw_alignment.fasta \
        -out_AA raw_alignment_AA.fasta > /dev/null 2>&1
fi

# 2. POST-ALIGNMENT FILTERING
# We use 'goalign' as it deals with alignments efficiently
echo "[2/4] Filtering Alignment..."

# Step 2a: Clean SITES (Columns)
# Remove columns where >50% of sequences have a gap. 
# This removes regions that are likely insertion artifacts in a few sequences.
goalign clean sites -c "$MAX_GAP_SITE" -i raw_alignment.fasta -o filtered_sites.fasta

# Step 2b: Clean SEQUENCES (Rows)
# Remove sequences that are >50% gaps (after site cleaning).
goalign clean seqs -c "$MAX_GAP_SEQ" -i filtered_sites.fasta -o "${OUTPUT_PREFIX}.fasta"

# 3. FINAL STATS
echo "[3/4] Generating Final Report..."
echo "--- Raw Alignment Stats ---"
seqkit stats raw_alignment.fasta
echo "--- Filtered Alignment Stats ---"
seqkit stats "${OUTPUT_PREFIX}.fasta"

# Optional: Convert to other formats if needed (e.g., Phylip for PAML)
# goalign reformat phylip -i "${OUTPUT_PREFIX}.fasta" -o "${OUTPUT_PREFIX}.phy"

echo "--- DONE. Final file: ${OUTPUT_PREFIX}.fasta ---"