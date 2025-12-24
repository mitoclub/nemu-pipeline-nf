#!/bin/bash

# =================CONFIGURATION=================
INPUT_FILE="$1"          
GENCODE="${2:-1}"        
OUTPUT_PREFIX="final_output"
THREADS=4
MACSE_JAR="/opt/macse_v2.07.jar"

# THRESHOLDS
LARGE_DATA_CUTOFF=500    
MIN_SEQ_LEN=100          
MAX_GAP_SEQ=0.50         
MAX_GAP_SITE=0.50        
# ===============================================

if [[ -z "$INPUT_FILE" ]]; then
    echo "Usage: ./align_and_filter.sh <input.fasta> [gencode]"
    exit 1
fi

echo "--- STARTING PIPELINE ---"

# [0] PRE-FILTERING
seqkit seq -m $MIN_SEQ_LEN -g "$INPUT_FILE" > input_clean.fasta
SEQ_COUNT=$(grep -c "^>" input_clean.fasta)

# [1] ALIGNMENT (Big Data vs Pure MACSE)
if [ "$SEQ_COUNT" -gt "$LARGE_DATA_CUTOFF" ]; then
    echo "[1/5] Strategy: BIG DATA (> $LARGE_DATA_CUTOFF sequences)"
    
    # 1a. Trim & Repair
    java -jar "$MACSE_JAR" -prog trimNonHomologousFragments -seq input_clean.fasta -gc_def "$GENCODE" -out_NT 1_trimmed.fasta > /dev/null 2>&1
    java -jar "$MACSE_JAR" -prog alignSequences -seq 1_trimmed.fasta -gc_def "$GENCODE" -out_NT 2_repaired.fasta -max_refine_iter 0 -local_realign_init 0 > /dev/null 2>&1
    sed 's/-//g' 2_repaired.fasta > 3_ungapped.fasta

    # 1b. Translate & Align
    java -jar "$MACSE_JAR" -prog translateNT2AA -seq 3_ungapped.fasta -gc_def "$GENCODE" -out_AA 4_protein.faa > /dev/null 2>&1
    mafft --thread "$THREADS" --auto --quiet 4_protein.faa > 5_aligned_protein.faa

    # 1c. Back-Translate
    java -jar "$MACSE_JAR" -prog reportGapsAA2NT -align_AA 5_aligned_protein.faa -seq 3_ungapped.fasta -gc_def "$GENCODE" -out_NT raw_alignment.fasta > /dev/null 2>&1

    rm 1_trimmed.fasta 2_repaired.fasta 3_ungapped.fasta 4_protein.faa 5_aligned_protein.faa
else
    echo "[1/5] Strategy: PURE MACSE (< $LARGE_DATA_CUTOFF sequences)"
    java -jar "$MACSE_JAR" -prog alignSequences -seq input_clean.fasta -gc_def "$GENCODE" -out_NT raw_alignment.fasta -out_AA raw_alignment_AA.fasta > /dev/null 2>&1
fi

# [2] SANITIZING (Export Alignment) - NEW STEP
# Replaces '!' and internal stops with 'NNN' or '---' to fix biological validity
echo "[2/5] Sanitizing Alignment (Masking stops/frameshifts)..."

java -jar "$MACSE_JAR" -prog exportAlignment \
    -align raw_alignment.fasta \
    -gc_def "$GENCODE" \
    -ambi_OFF \
    -codonForInternalStop "NNN" \
    -codonForFinalStop "---" \
    -codonForInternalFS "NNN" \
    -codonForExternalFS "---" \
    -out_NT sanitized_alignment.fasta \
    -out_stat_per_seq macse_stat_per_seq.csv \
    -out_stat_per_site macse_stat_per_site.csv > /dev/null 2>&1

# [3] FILTERING (GoAlign)
echo "[3/5] Filtering Alignment..."

# 3a. Clean Sites (Columns)
goalign clean sites -c "$MAX_GAP_SITE" -i sanitized_alignment.fasta -o filtered_sites.fasta

# 3b. Clean Sequences (Rows)
goalign clean seqs -c "$MAX_GAP_SEQ" -i filtered_sites.fasta -o "${OUTPUT_PREFIX}.fasta"

# [4] FINAL STATS
echo "[4/5] Generating Final Report..."
echo "--- Final Filtered Stats ---"
seqkit stats "${OUTPUT_PREFIX}.fasta"

# [5] OPTIONAL: N-Content Warning
# If a sequence is mostly 'NNN' (from your sanitizer), you might want to know.
echo "--- Ambiguity Check (Sequences with >20% Ns) ---"
seqkit fx2tab --name --gc --avg-qual "${OUTPUT_PREFIX}.fasta" | awk '$4 > 20 {print $1 " has high N content"}' 

echo "--- DONE. Final file: ${OUTPUT_PREFIX}.fasta ---"