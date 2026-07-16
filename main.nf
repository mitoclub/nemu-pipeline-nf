#!/usr/bin/env nextflow

/*
 * Cite as:
 * NeMu: a comprehensive pipeline for accurate reconstruction of neutral mutation spectra from evolutionary data
 * https://doi.org/10.1093/nar/gkae438
 * Authors: Bogdan Efimenko, Konstantin Popadin, Konstantin Gunbin
 */

/* Requirements:
 * Nextflow 25.10, seqkit, taxonkit, BLAST+, Python 3.14, mafft, macse, goalign, iqtree2, newick_utils
 * pymutspec 0.0.15 (Python lib)
 */

// --- Global Parameters ---

// Inputs/Outputs
params.input            = ""

// Databases & Tools
params.db               = ""
params.taxdump          = ""

// Pipeline Logic
        // TODO rename nucleotide_... to short versions (nuc_cds and nuc_any)
params.inputType        = "cds"                         // protein, cds, noncds
params.speciesName      = ""                            // Override species name
params.gencode          = 1
params.maxTargetSeqs    = 2000
params.minSeqs          = 4                             // Min sequences to proceed
params.threads          = 1
params.saveIntermeds    = false                         // Save intermediate files TODO  
params.help             = false                          // Show help message
params.verbose          = false

// MSA & Tree
params.outgroupId       = "OUTGRP"                      // Manually specify outgroup sequence ID for 'nucleotide' input
params.aligned          = false                         // 'nucleotide' input can be pre-aligned
params.msaMode          = "auto"                        // auto, macse (accurate codon alignment), mafft_macse (fast codon alignment), mafft
params.treefile         = ""                            // Provide pre-computed tree path to skip tree building
params.model            = "GTR+FO+G4+I"                 // IQ-TREE Model
params.modelAsr         = "GTR+FO+G4+I"                 // ASR Model
params.runTreeShrink    = true                          // Run TreeShrink to prune long branches

// Mutation Extraction
params.consCatCutoff    = 0                             // 0 = no cutoff // TODO pass list of categories instead of single value
params.probaArg         = true                          // Use probabilities
params.uncertaintyCoef  = true                          // Use phylogeny uncertainty coefficient TODO improve implementation

// Spectra Calculation
params.plot             = false                         // Generate plots

// Subsets of mutations to derive spectra for
params.internal         = false                         // TODO rename and terminal too
params.terminal         = false
params.branchSpectra    = false

// Type of spectra that will be aggregated to main output table 'spectra_total.tsv'. Can be 'syn', 'syn4f', 'nonsyn', 'all'
params.spectraType      = "syn"  // TODO implement
params.calc192          = true   // TODO implement

process PREPARE_TAXONOMY {
    input:
    path species_list
    path taxdump_dir

    output:
    path "parsed_taxonomy.txt"

    script:
    """
    export TAXONKIT_DB=${taxdump_dir}

    # derive the taxids for input species names, but also keep numeric taxids if they were provided directly
    taxonkit name2taxid $species_list | awk -F'\t' -v OFS="\t" '\$1 ~ /^[0-9]+\$/ { \$2 = \$1 } { print }' > species_taxids.txt

    if [ ! -s species_taxids.txt ]; then
        : > lineages.txt
        echo '{}' > taxonlist.json
        : > parsed_taxonomy.txt
        exit 0
    fi

    taxonkit lineage --taxid-field 2 species_taxids.txt | taxonkit reformat -t -f "{f},{s}" --lineage-field 3 > lineages.txt

    # Prepare a single id list for taxonkit list --ids.
    tx_list_to_parse=\$(cut -f5 lineages.txt | tr ',' '\n' | awk '/^[0-9]+\$/' | sort -nu | paste -sd ',' -)
    if [ -z "\$tx_list_to_parse" ]; then
        tx_list_to_parse=\$(cut -f2 species_taxids.txt | paste -sd ',' -)
    fi

    if [ -n "\$tx_list_to_parse" ]; then
        taxonkit list --ids "\$tx_list_to_parse" --json > taxonlist.json
    else
        echo '{}' > taxonlist.json
    fi

    # extract paths to leaf nodes (species) for lineage parsing in the next step
    cat taxonlist.json | jq -r 'paths(objects | select(length == 0)) | join(",")' > paths.csv

    cut -f5 lineages.txt > fam_sp_taxids.csv

    : > descendants.txt

    # iterate over family and species taxids and extract their lineages from taxonlist.json
    mkdir -p sp_lineage fam_lineage fam_lineage_excl_sp
    while IFS=, read -r fam_taxid sp_taxid; do
        echo "Processing \$fam_taxid and \$sp_taxid"
        echo \$sp_taxid > sp_lineage/\${sp_taxid}.txt
        cat taxonlist.json | jq --arg tax \$sp_taxid -r '.[\$tax]  | paths(objects | select(length == 0)) | join("\n")' | sort -nu >> sp_lineage/\${sp_taxid}.txt
        cat taxonlist.json | jq --arg tax \$fam_taxid -r '.[\$tax] | paths(objects | select(length == 0)) | join("\n")' | sort -nu > fam_lineage/\${fam_taxid}.txt
        grep -vf sp_lineage/\${sp_taxid}.txt fam_lineage/\${fam_taxid}.txt > fam_lineage_excl_sp/\${fam_taxid}.txt

        sp_lineage_lst=\$(paste -sd ',' sp_lineage/\${sp_taxid}.txt)
        fam_lineage_lst=\$(paste -sd ',' fam_lineage_excl_sp/\${fam_taxid}.txt)

        paste <(echo "\$sp_lineage_lst") <(echo "\$fam_lineage_lst") >> descendants.txt

    done < fam_sp_taxids.csv
    
    # combine lineages and descendants into final taxonomy file if their line numbers match
    if [ \$(wc -l < lineages.txt) -ne \$(wc -l < descendants.txt) ]; then
        echo "WARNING: Line count mismatch between lineages and descendants. Check the intermediate files for details."
        exit 1
    else
        paste <(cat lineages.txt) <(cat descendants.txt) >> parsed_taxonomy.txt
    fi
    """
}

process TBLASTN {
    tag "$id"
    cpus params.threads
    errorStrategy 'ignore'
    
    input:
    tuple val(id), val(species_name), val(sequence)
    val db_path
    path parsed_taxonomy
    val gencode
    val max_target_seqs

    output:
    tuple val(id), path("blast_species.tsv"), path("blast_outgroup.tsv")

    script:
    outfmt="6 saccver pident length qlen gapopen sstart send evalue bitscore sframe"
    """
    # Save query
    echo ">${id}" > query.fa
    echo "${sequence}" >> query.fa

    # extract species and family taxids for the given species name from parsed_taxonomy
    awk -F'\t' -v OFS="\t" -v name="$species_name" '
        \$1 == name {print}
    ' $parsed_taxonomy > matched_species.txt

    if [ ! -s matched_species.txt ]; then
        echo "No matching species found in taxonomy for ${species_name}. Skipping BLAST."
        touch blast_species.tsv blast_outgroup.tsv
        exit 0
    fi

    cut -f6 matched_species.txt | tr "," "\n" | sort -n > species_taxids.txt
    cut -f7 matched_species.txt | tr "," "\n" > family_taxids.txt
    
    # 1. Species BLAST
    if [ -s species_taxids.txt ]; then
        echo "Running Species BLAST..."
        tblastn -query query.fa -db ${db_path} -db_gencode ${gencode} \
            -max_target_seqs ${max_target_seqs} \
            -evalue 0.0001 -num_threads ${task.cpus} \
            -taxidlist species_taxids.txt -no_taxid_expansion \
            -outfmt "$outfmt" \
            -out blast_species.tsv || touch blast_species.tsv
    else
        echo "No species taxids found. Skipping Species BLAST."
        touch blast_species.tsv
    fi

    # 2. Outgroup BLAST (Relatives)
    if [ -s family_taxids.txt ]; then
        echo "Running Outgroup BLAST..."
        # Only need top 10 hits to find a good outgroup
        tblastn -query query.fa -db ${db_path} -db_gencode ${gencode} \
            -max_target_seqs 10 \
            -evalue 0.0001 -num_threads ${task.cpus} \
            -taxidlist family_taxids.txt -no_taxid_expansion \
            -outfmt "$outfmt" \
            -out blast_outgroup.tsv || touch blast_outgroup.tsv
    else
        echo "No relative taxids found. Skipping Outgroup BLAST."
        touch blast_outgroup.tsv
    fi
    """
}

process FILTER_AND_EXPORT {
    tag "$id"
    cpus params.threads
    errorStrategy 'ignore'

    input:
    tuple val(id), path(records_species), path(records_relatives)
    val db_path

    output:
    tuple val(id), path("sampled_sequences.fasta"), env("OUTGRP_ID"), emit: seqs
    path "blast_filtering_log.txt"

    script:
    """
    python3 -c "
print('INFO: Filtering BLAST hits')
hits = []
outgroup_hits = []

# Thresholds
MIN_IDENT_SPECIES = 80.0
MIN_COV_SPECIES = 0.5
MIN_IDENT_OUTGROUP = 70.0 
MIN_COV_OUTGROUP = 0.5

def parse_blast(filename, hits_list, min_ident, min_cov):
    try:
        with open(filename) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 10: continue
                sacc, pident = parts[0], float(parts[1])
                length, qlen = float(parts[2]), float(parts[3])
                bitscore, sframe = float(parts[8]), int(parts[9])
                sstart, send = parts[5], parts[6]
                coverage = length / qlen
                
                if pident >= min_ident and coverage >= min_cov:
                    hits_list.append({
                        'sacc': sacc, 'sstart': sstart, 'send': send, 
                        'sframe': sframe, 'bitscore': bitscore, 'pident': pident
                    })
    except FileNotFoundError:
        pass

parse_blast('$records_species',  hits,  MIN_IDENT_SPECIES,  MIN_COV_SPECIES)
parse_blast('$records_relatives', outgroup_hits, MIN_IDENT_OUTGROUP, MIN_COV_OUTGROUP)

if outgroup_hits:
    outgroup_hits.sort(key=lambda x: x['bitscore'], reverse=True)
    best_out = outgroup_hits[0]
    print(f\\"Selected Outgroup: {best_out['sacc']}, Ident: {best_out['pident']}%, Score: {best_out['bitscore']}\\")
    hits.append(best_out)    
else:
    print('WARNING: No valid outgroup found.')

print(f'Found {len(hits)} valid species hits.')
final_coords = []
for h in hits:
    # Logic: ID start-end strand
    strand = 'minus' if h['sframe'] < 0 else 'plus'
    x, y = (h['sstart'], h['send']) if strand == 'plus' else (h['send'], h['sstart'])
    final_coords.append(f\\"{h['sacc']} {x}-{y} {strand}\\")

with open('extract_coords.txt', 'w') as f:
    for line in final_coords:
        f.write(line + '\\n')
    " >> blast_filtering_log.txt
    # END OF PYTHON CODE

    OUTGRP_ID=\$(grep -oP "Selected Outgroup: \\K[^,]+" blast_filtering_log.txt)

    # 4. Extract
    if [ -s extract_coords.txt ]; then
        blastdbcmd -db ${db_path} -entry_batch extract_coords.txt -outfmt %f -out sampled_sequences.fasta
    else
        touch sampled_sequences.fasta
    fi
    """
}

process ENCODE_AND_RMDUP {
    tag "$id"
    errorStrategy 'ignore'

    input:
    tuple val(id), path(sequences), val(OUTGRP_ID)

    output:
    tuple val(id), path("seqs_unique.fasta"), path("encoded_headers.txt"), env('NUM_SEQS')

    script:
    """
    # TODO error for nucl input: OUTGRP_ID not specified, it in description
    seqkit replace -p .+ -r "seq_{nr}" -w 0 < $sequences > encoded_raw.fasta
    if [ -z "${OUTGRP_ID}" ]; then
        mv encoded_raw.fasta encoded.fasta
    else
        # this logic works only for protein input
        # TODO find line with outgrp sign and replace its basid id (seq_i) with OUTGRP
        # !!!! align old and encoded headers/ids and replace correctly
        outgrp_id=\$(seqkit seq -i -n < ./encoded_raw.fasta | tail -1)
        seqkit replace -p \${outgrp_id} -r "OUTGRP" -w 0 < encoded_raw.fasta > encoded.fasta
    fi

    # Save mapping
    codes=\$(seqkit seq -ni < ./encoded.fasta)
    original_names=\$(seqkit seq -n < ./${sequences})
    paste <(echo "\$codes") <(echo "\$original_names") > encoded_headers.txt
    
    # Remove duplicates
    seqkit rmdup -D duplicated.txt -s -w 0 < encoded.fasta > seqs_unique.fasta
    NUM_SEQS=\$(grep -c '>' ./seqs_unique.fasta)
    """
}

process MSA {
    tag "$id"
    cpus params.threads
    errorStrategy 'ignore'

    input:
    tuple val(id), path(sequences), path(encoded_headers), val(num_seqs)
    val gencode
    val msa_mode

    output:
    tuple val(id), path("msa.fasta"), path("alignment.log"), env('NUM_SEQS_FLT')

    script:
    // THRESHOLDS
    LARGE_DATA_CUTOFF=50    // Switch to Big Data workflow if seqs > this
    MIN_SEQ_LEN=100          // Pre-filter: remove sequences shorter than 100bp
    MAX_GAP_SEQ=0.50         // Post-filter: Remove SEQS with >50% gaps
    MAX_GAP_SITE=0.50        // Post-filter: Remove SITES (columns) with >50% gaps

    """
    echo "--- STARTING MSA ---" > alignment.log
    echo "[0/5] Pre-filtering short sequences (<${MIN_SEQ_LEN}bp)..." >> alignment.log
    seqkit seq -m $MIN_SEQ_LEN -g "$sequences" > input_clean.fasta
    SEQ_COUNT=\$(grep -c "^>" input_clean.fasta)
    echo "Sequences remaining: \$SEQ_COUNT" >> alignment.log

    # MSA Mode Selection
    if [ $msa_mode = "auto" ]; then
        if [ "\$SEQ_COUNT" -gt "$LARGE_DATA_CUTOFF" ]; then
            msa_mode_sh="mafft_macse"
        else
            msa_mode_sh="macse"
        fi
    else
        msa_mode_sh="$msa_mode"
    fi

    # ALIGNMENT
    if [ \$msa_mode_sh = "mafft_macse" ]; then
        # Big Data Strategy (Trim -> MacseRepair -> Mafft -> MacseBackTrans)
        macse -prog trimNonHomologousFragments \
            -seq input_clean.fasta -gc_def "$gencode" \
            -out_NT 1_trimmed.fasta

        macse -prog alignSequences \
            -seq 1_trimmed.fasta -gc_def "$gencode" \
            -out_NT 2_repaired.fasta \
            -max_refine_iter 0 -local_realign_init 0

        sed 's/-//g' 2_repaired.fasta > 3_ungapped.fasta

        macse -prog translateNT2AA \
            -seq 3_ungapped.fasta -gc_def "$gencode" \
            -out_AA 4_protein.faa

        mafft --thread "${task.cpus}" --auto --quiet 4_protein.faa > 5_aligned_protein.faa

        macse -prog reportGapsAA2NT \
            -align_AA 5_aligned_protein.faa \
            -seq 3_ungapped.fasta -out_NT raw_alignment.fasta

        # Cleanup intermediate files # TODO enable if needed
        # rm 1_trimmed.fasta 2_repaired.fasta 3_ungapped.fasta 4_protein.faa 5_aligned_protein.faa

    elif [ \$msa_mode_sh = "macse" ]; then
        # Pure MACSE
        macse -prog alignSequences \
            -seq input_clean.fasta -gc_def "$gencode" \
            -out_NT raw_alignment.fasta \
            -out_AA raw_alignment_AA.fasta
    
    elif [ \$msa_mode_sh = "mafft" ]; then
        mafft --thread "${task.cpus}" --auto --quiet input_clean.fasta > raw_alignment.fasta
    fi

    # SANITIZING
    if [ \$msa_mode_sh != "mafft" ]; then
        macse -prog exportAlignment \
            -align raw_alignment.fasta \
            -gc_def "$gencode" \
            -ambi_OFF \
            -codonForInternalStop "NNN" -codonForFinalStop "---" \
            -codonForInternalFS "NNN" -codonForExternalFS "---" \
            -out_NT sanitized_alignment.fasta \
            -out_stat_per_seq macse_stat_per_seq.csv \
            -out_stat_per_site macse_stat_per_site.csv
    else
        mv raw_alignment.fasta sanitized_alignment.fasta
    fi

    # POST-ALIGNMENT FILTERING
    # Remove columns where >50% of sequences have a gap. 
    # Remove sequences that are >50% gaps (after site cleaning).
    goalign clean sites -c "$MAX_GAP_SITE" -i sanitized_alignment.fasta -o filtered_sites.fasta
    goalign clean seqs -c "$MAX_GAP_SEQ" -i filtered_sites.fasta -o "filtered_seqs.fasta"
    seqkit rmdup -s < filtered_seqs.fasta > msa.fasta

    echo "--- Alignments Stats ---" >> alignment.log
    seqkit stats raw_alignment.fasta msa.fasta >> alignment.log
    
    # N-Content Warning
    seqkit fx2tab --name --gc --avg-qual "msa.fasta" | \
        awk '\$4 > 20 {print \$1 " has high N content"}' >> alignment.log
    NUM_SEQS_FLT=\$(grep -c "^>" msa.fasta)
    """
}


process BUILD_TREE {
    tag "$id"
    cpus params.threads
    errorStrategy 'ignore'
    
    input:
    tuple val(id), path(sequences)
    val model
    val run_treeshrink

    output:
    tuple val(id), path("msa_filtered.fasta"), path("tree_rerooted.nwk")

    script:
    QUANTILE=0.1
    """
    echo "Building tree de novo..."
    iqtree2 -s $sequences -m $model -nt $task.cpus --prefix ml

    # Set OUTGRP_PARAM to "-x OUTGRP" if outgroup is present, else empty string
    if grep -q '^>OUTGRP' "$sequences"; then
        OUTGRP_PARAM="-x OUTGRP"
    else
        OUTGRP_PARAM=""
    fi

    nseq=\$(grep -c '>' $sequences)
    if [ $run_treeshrink = true ] && [ \$nseq -gt 10 ]; then
        run_treeshrink.py -t ml.treefile -O treeshrink -o . -q $QUANTILE \$OUTGRP_PARAM
        mv treeshrink.treefile treeshrink.nwk
    else
        mv ml.treefile treeshrink.nwk
    fi

    # Determine if OUTGRP is present in the resulting tree
    if nw_labels treeshrink.nwk | grep -qx 'OUTGRP'; then
        # Check outgroup "quality"
        nw_distance -m p -s f -n treeshrink.nwk > branches.txt
        LC_ALL=C sort -grk 2 branches.txt > branches.txt.sorted
        
        # Prune bad outgroup if needed (simple heuristic: if branch to OUTGRP has not so large length)
        # TODO get 10% of branches
        head -n 5 branches.txt.sorted > branches.txt.top5
        if grep -q OUTGRP branches.txt.top5; then
            nw_reroot -l treeshrink.nwk OUTGRP > tree_rerooted.nwk
        else
            nw_prune treeshrink.nwk OUTGRP | nw_reroot - > tree_rerooted.nwk
        fi
    else
        # No explicit outgroup present: root on the longest branch
        nw_reroot -l treeshrink.nwk > tree_rerooted.nwk
    fi

    # drop sequences not present in the tree
    nw_labels -I tree_rerooted.nwk  > leaves.txt
    seqkit grep -f leaves.txt -w 0 $sequences > msa_filtered.fasta
    seqkit stats $sequences msa_filtered.fasta # sanity check
    """
}

process INCLUDE_USER_TREE {
    tag "$id"
    errorStrategy 'ignore'

    input:
    tuple val(id), path(sequences)
    path treefile 

    output:
    tuple val(id), path("msa_filtered.fasta"), path(treefile)

    script:
    """
    # drop sequences not present in the tree
    nw_labels -I "$treefile"  > leaves.txt
    seqkit grep -f leaves.txt -w 0 $sequences > msa_filtered.fasta
    """
}

process ASR {
    tag "$id"
    cpus params.threads
    errorStrategy 'ignore'

    input:
    tuple val(id), path(sequences), path(tree)
    val model

    output:
    tuple val(id), path(sequences), path("tree.nwk"), path("iqtree_anc.state"), path("rates.tsv")

    script:
    """
    iqtree2 -te $tree -s $sequences -m $model -asr -nt $task.cpus --prefix asr --rate
    mv asr.rate rates.tsv

    if grep -q OUTGRP asr.treefile; then
        nw_reroot -l asr.treefile OUTGRP | sed 's/;/ROOT;/' > tree.nwk
    else
        nw_reroot asr.treefile | sed 's/;/ROOT;/' > tree.nwk
    fi

    iqtree_states_add_part.py asr.state iqtree_anc.state
    """
}

process DRAW_TREE {
    tag "$id"

    input:
    tuple val(id), path(tree), path(nodes_mapping)

    output:
    tuple val(id), path("tree.svg"), path("tree.png")

    script:
    """
    cut -f2 "$nodes_mapping" | sed 's/[^A-Za-z0-9_]/_/g' | cut -c 1-20 > cleaned_headers.txt
    paste <(cut -f1 "$nodes_mapping") cleaned_headers.txt > nodes_mapping_cleaned.txt
    
    nw_rename $tree nodes_mapping_cleaned.txt > tree_renamed.nwk
    nw_display -s -b 'visibility:hidden' -i 'visibility:hidden' tree_renamed.nwk > tree.svg
    
    # Convert SVG -> PNG
    rsvg-convert tree.svg -o tree.png || python3 -c "
try:
    from cairosvg import svg2png
    svg2png(url='tree.svg', write_to='tree.png')
except Exception as e:
    print('ERROR: No SVG->PNG converter found (rsvg-convert or cairosvg).')
    "
    """
}

process MUT_EXTRACTION {
    tag "$id"
    cpus params.threads
    errorStrategy 'ignore'

    input:
    tuple val(id), path(sequences), path(tree), path(internal_states), path(rates)
    val gencode
    val proba_arg
    val uncertainty_coef
    val cons_cat_cutoff 

    output:
    tuple val(id), path("observed_mutations.tsv"), path("expected_freqs.tsv"), emit: mutations
    // tuple val(id), path(sequences), path(tree), path("mut_extraction.log")
    tuple val(id), path("mut_extraction.log"), emit: logs
    // path "expected_mutations.tsv.gz"

    script:
    """
    # --threads ${task.cpus}
    ARGS="--gencode $gencode --no-mutspec --outdir mout --syn --syn4f --nonsyn --save-exp-muts"
    if [ $cons_cat_cutoff -gt 0 ]; then
        ARGS="\$ARGS --rates $rates --cat-cutoff $cons_cat_cutoff"
    fi
    if [ $proba_arg = "true" ]; then
        ARGS="\$ARGS --proba --pcutoff 0.3"
    fi
    if [ $uncertainty_coef = "true" ]; then
        ARGS="\$ARGS --phylocoef"
    else
        ARGS="\$ARGS --no-phylocoef"
    fi

    alignment2iqtree_states.py $sequences terminal_states.state

    collect_mutations.py --tree $tree --states terminal_states.state --states $internal_states \
        \$ARGS 

    mv mout/* .
    mv mutations.tsv observed_mutations.tsv
    mv run.log mut_extraction.log
    gzip expected_mutations.tsv
    """
}

process DERIVE_SPECTRA {
    tag "$id"
    errorStrategy 'ignore'

    input:
    tuple val(id), path(obs_muts), path(exp_freqs)
    val plot
    val internal
    val terminal
    val branch_spectra

    output:
    path "ms12syn_labeled.txt", emit: syn_spectrum
    tuple val(id), path("*.tsv"), emit: spectra_data
    tuple val(id), path("*.png"), optional: true, emit: spectra_plots

    script:
    """
    nmuts=`cat $obs_muts | wc -l`
    if [ \$nmuts -lt 2 ]; then
        echo "ERROR: There are no reconstructed mutations." >&2
        exit 1 # TODO handle this better check line num in the nextflow level
    fi

    ARGS="--exclude OUTGRP,ROOT --mnum192 16 --proba_cutoff 0.3 --syn --syn4f --all --nonsyn"
    if [ $plot = true ]; then
        ARGS="\$ARGS --plot -x png"
    fi

    # TODO replace mean_expected_mutations.tsv with exp_muts if needed
    # TODO replace by pure python code

    # TODO fix error "After filtration observed 0 mutations" 36__K10914
    
    # Main Calculation
    calculate_mutspec.py -b $obs_muts -e $exp_freqs -o . \$ARGS
    
    if [ ! -f ms12syn.tsv ]; then
        # TODO improve filtration quality in the script
        touch ms12syn_labeled.txt
        exit 0
    fi

    # Internal
    if [ "$internal" = "true" ]; then
        calculate_mutspec.py -b $obs_muts -e $exp_freqs -o . \$ARGS --subset internal
        # Cleanup
        if [ -f mean_expexted_mutations_internal.tsv ]; then rm mean_expexted_mutations_internal.tsv; fi
    fi
    
    # Terminal
    if [ "$terminal" = "true" ]; then
        calculate_mutspec.py -b $obs_muts -e $exp_freqs -o . \$ARGS --subset terminal
        if [ -f mean_expexted_mutations_terminal.tsv ]; then rm mean_expexted_mutations_terminal.tsv; fi
    fi
    
    # Branch Spectra
    if [ "$branch_spectra" = "true" ]; then
        calculate_mutspec.py -b $obs_muts -e $exp_freqs -o . \$ARGS --branches
    fi
     
    # add column $id to output tsv file
    python3 -c "
import pandas as pd
df = pd.read_csv('ms12syn.tsv', sep='\\t')
df['QueryId'] = '${id}'
df.to_csv('ms12syn_labeled.txt', sep='\\t', index=False)
    "
    """
}

// process AMINO_ACID_STUFF TODO

process CHECK_INPUT_TYPE {
    input:
    path fasta

    output:
    tuple path(fasta), env("TYPE")

    script:
    """
    TYPE=\$(seqkit stats $fasta -T | tail -1 | cut -f3)

    # TODO uppercase with seqkit
    # seqkit seq 
    """
}

process WRITE_README {
    output:
    path("readme.txt")

    script:
"""
cat > readme.txt <<- EOM
Output structure:

TODO update after all changes

outputDir
├── parsed_id_name
│   ├── encoded_headers.txt             # Mapping of encoded headers to original headers
│   ├── expected_freqs.tsv              # 
│   ├── images
│   │   ├── ms12all.png
│   │   ├── ms12ff.png
│   │   ├── ms12nonsyn.png
│   │   ├── ms12syn.png
│   │   ├── ms192all.png
│   │   ├── ms192ff.png
│   │   ├── ms192syn.png
│   │   ├── tree.png
│   │   └── tree.svg
│   ├── iqtree_anc.state
│   ├── mean_expexted_mutations.tsv
│   ├── ms12all.tsv
│   ├── ms12ff.tsv
│   ├── ms12nonsyn.tsv
│   ├── ms12syn.tsv
│   ├── ms192all.tsv
│   ├── ms192ff.tsv
│   ├── ms192syn.tsv
│   ├── msa_filtered.fasta
│   ├── mut_extraction.log
│   ├── observed_mutations.tsv
│   ├── rates.tsv
│   ├── seqs_unique.fasta
│   └── tree.nwk
├── readme.txt
└── spectra_total.tsv

EOM
"""
}

// Helper: Check dependencies
boolean commandExists(String command) {
    def proc = ["bash", "-c", "command -v ${command}"].execute()
    proc.waitFor()
    return proc.exitValue() == 0
}

def printHelpMessage(String version, params) {
    println """
N E M U   P I P E L I N E  ${version}
================================
A comprehensive pipeline for accurate reconstruction of neutral mutation spectra from evolutionary data.
https://doi.org/10.1093/nar/gkae438

Usage: nextflow run main.nf --input <input_fasta> [options]

TODO update to latest

Main options:
    --input FILE            Input FASTA file (required)
                            If input type is protein, a multi-FASTA with one or several sequences 
                            required (header format: ">ID [Species name]"). If input type 
                            is nucleotide, one or several fasta files with orthologous 
                            sequences (including outgroup) required
    --input_type STRING     Type of input sequences: 
                            protein, cds, noncds (default: cds)
    --gencode NUM           Genetic code table (default: 1)
                            Used for codon-aware alignment and annotation of mutations

Nextflow options:
    -with-report FILE       Generate execution report
    -with-trace FILE        Generate execution trace
    -with-timeline FILE     Generate execution timeline
    -output-dir DIR         Specify output directory (default: ./results)

Common options:
    --threads NUM           Number of threads to use (default: 1) TODO update description
    --save-intermeds BOOL   Save intermediate files TODO implement
    --help                  Show this help message and exit

Options for nucleotide input:
    --outgroup-id STRING    Outgroup sequence ID (default: ${params.outgroupId})
                            Specify outgroup sequence ID in the alignment for rooting the tree
    --aligned BOOL          Input sequences are pre-aligned (default: ${params.aligned})

Options for protein input:
    --db PATH               BLAST database path (required for protein input)
    --taxdump DIR           Taxdump directory path (required for protein input)
                            TODO integrate to the container
    --species-name STRING   Override species name. Useful when you work with proteins 
                            from single species
    --max-target-seqs NUM   Max target sequences for BLAST (default: ${params.maxTargetSeqs})
                            tblastn will collect no more than this number of sequences
    
Options for MSA & Phylogeny:
    --msa-mode STRING       MSA mode: auto, macse, mafft_macse, mafft (default: ${params.msaMode})
    --min-seqs NUM          Minimum number of sequences to proceed phylogenetic inference (default: ${params.minSeqs})
    --treefile FILE         Input tree file (optional; default: build tree de novo)
    --model STRING          IQ-TREE substitution model (default: ${params.model})
    --model-asr STRING      ASR substitution model (default: ${params.modelAsr})
    --run-treeshrink BOOL   Run TreeShrink to prune long branches (default: ${params.runTreeShrink})

Options for Mutation Extraction:
    --cons-cat-cutoff NUM   Conservation category cutoff for mutation extraction (default: ${params.consCatCutoff})
                            0 = no cutoff; only mutations in sites with rate category 
                            less than or equal to this value will be used
    --proba-arg BOOL        Use probabilities in mutation extraction (default: ${params.probaArg})
    --uncertainty-coef BOOL Use phylogeny uncertainty coefficient in mutation extraction 
                            (default: ${params.uncertaintyCoef})

Options for Mutation Spectra Derivation:
    --plot BOOL             Generate barcharts of mutation spectra (default: ${params.plot})
    --internal BOOL         Derive spectra for internal branches (default: ${params.internal})
    --terminal BOOL         Derive spectra for terminal branches (default: ${params.terminal})
    --branch-spectra BOOL   Derive spectra for individual branches (default: ${params.branchSpectra})
""".stripIndent()
}


workflow nemuCore {
    take:
    nucl_multi_fasta
    gencode
    min_seqs
    aligned
    msa_mode
    model
    model_asr
    treefile
    run_treeshrink
    proba_arg
    uncertainty_coef
    cons_cat_cutoff
    plot
    internal
    terminal
    branch_spectra

    main:
    ENCODE_AND_RMDUP(nucl_multi_fasta)

    // Filter Low Count
    seq_num_verified_ch = ENCODE_AND_RMDUP.out.filter { id, _seq, _enc_head, num_seqs ->
        if (num_seqs.toInteger() > min_seqs) return true
        log.warn "SKIPPING ${id}: Count ${num_seqs} < ${min_seqs}"
        return false
    }

    if (aligned == true) {
        msa_num_verified_ch = seq_num_verified_ch.map { id, fasta, _enc_head, _num_seqs ->
            [id, fasta]
        }
    } else {
        MSA(seq_num_verified_ch, gencode, msa_mode)

        msa_num_verified_ch = MSA.out.filter { id, _seq, _msa_log, num_seqs ->
            if (num_seqs.toInteger() > min_seqs) return true
            log.warn "SKIPPING ${id}: Count after MSA ${num_seqs} < ${min_seqs}"
            return false
        }.map { id, fasta, _msa_log, _num_seqs -> [id, fasta] }
    }

    // Build or use user-provided tree
    if (treefile && treefile != "") {
        if (file(treefile).exists()) {
            log.info "Using user-provided tree file: ${treefile}"
            tree_ch = INCLUDE_USER_TREE(msa_num_verified_ch, treefile)
        } else {
            log.warn "User-provided tree file '${treefile}' does not exist. Ignoring and building tree de novo."
            tree_ch = BUILD_TREE(msa_num_verified_ch, model, run_treeshrink)
        }
    } else {
        tree_ch = BUILD_TREE(msa_num_verified_ch, model, run_treeshrink)
    }

    ASR(tree_ch, model_asr)

    // make channel with tree from ASR and nodes mapping from ENCODE_AND_RMDUP for drawing
    just_tree_ch = ASR.out.map { id, _seq, tree, _anc_state, _rates -> [id, tree] }
                          .join(seq_num_verified_ch.map { id, _seq, enc_head, _num_seqs -> [id, enc_head] })
    DRAW_TREE(just_tree_ch)

    MUT_EXTRACTION(ASR.out, 
        gencode, proba_arg, uncertainty_coef,
        cons_cat_cutoff,
    )

    DERIVE_SPECTRA(MUT_EXTRACTION.out.mutations,
        plot, internal, terminal, branch_spectra
    )

    emit:
    syn_spectrum = DERIVE_SPECTRA.out.syn_spectrum
    encoded_headers = seq_num_verified_ch // encoded sequences after rmdp
    msa_tree = ASR.out
    tree_images = DRAW_TREE.out
    mutation_logs = MUT_EXTRACTION.out.logs
    mutation_data = MUT_EXTRACTION.out.mutations
    spectra_data = DERIVE_SPECTRA.out.spectra_data
    spectra_plots = DERIVE_SPECTRA.out.spectra_plots
}

workflow blastHead {
    take:
    input_fasta
    species_name
    taxdump
    db
    max_target_seqs
    gencode

    main:
    def seq_counter = 0
    raw_sequences = channel.fromPath(input_fasta)
        .splitFasta(record: [id: true, header: true, seqString: true])
        .filter { record ->
            def seq = record.seqString.toUpperCase()
            if (seq.length() < 30) {
                log.warn "SKIPPING ${record.id}: Sequence too short (${seq.length()})."
                return false
            }
            if (seq =~ /[EFILPQZ]/) return true
            def dna_count = seq.count('A') + seq.count('C') + seq.count('G') + seq.count('T') + seq.count('N')
            if ((dna_count / seq.length()) > 0.95) {
                log.warn "SKIPPING ${record.id}: Looks like DNA."
                return false
            }
            return true
        }
        .map { record ->
            def count = ++seq_counter
            def clean_original_id = record.id.split()[0].replaceAll(/[^a-zA-Z0-9\.]/, '_')
            // parse species name from header if not provided by user
            // species name patterns: ">ID gene name [species name]"
            def species = species_name
            if (!species_name) {
                def sp_match = record.header =~ /\[([^\]]+)\]/
                if (sp_match) {
                    species = sp_match[0][1]
                } else {
                    log.warn "SKIPPING ${record.id}: Species name not found in header. It must be provided in square brackets [Species name]."
                    species = "unknown_species"
                }
            }
            def species_cleaned = species.replaceAll(/[^a-zA-Z0-9\_\-]/, '_')
            def unique_id = "${count}__${clean_original_id}__${species_cleaned}"
            [unique_id, species, record.seqString]
        }.filter { _id, species, _seq ->
            return species != "unknown_species"
        }

    // show parsed sequences
    if (params.verbose) {
        raw_sequences.view { id, species, seq ->
            "Parsed sequence: ${id}, species: ${species}, length: ${seq.length()}"
        }
    }

    species_lst = raw_sequences.map { _id, species, _seq -> [species] }
        .unique().flatten().collectFile(name: 'species_lst.txt', newLine: true, sort: true)
    
    if (params.verbose) {
    species_lst.subscribe { file ->
            println "Parsed species names are saved to file: $file\n\n"
            // println "File content is:\n${file.text}"
        }
    }
            
    PREPARE_TAXONOMY(species_lst, taxdump)
    parsed_taxonomy = PREPARE_TAXONOMY.out.first()

    TBLASTN(raw_sequences, db, parsed_taxonomy, gencode, max_target_seqs)
    records = TBLASTN.out.filter { id, rec_sp, _rec_rel ->
        if (rec_sp.size() > 0) return true
        log.warn "SKIPPING ${id}: Species records not found."
        return false
    }
    FILTER_AND_EXPORT(records, db)
    
    // Filter missing sequences
    fasta = FILTER_AND_EXPORT.out.seqs.filter { id, fasta, outgrp_id ->
        if (outgrp_id == '') log.warn "Outgroup not found for ${id}. Continue anyway."
        if (fasta != null && fasta.size() > 0) return true
        log.warn "SKIPPING ${id}: No sequences found."
        return false
    }

    emit:
    nuc_multifasta = fasta
    taxonomy = parsed_taxonomy
}

workflow {
    main:
    NEMU_VERSION="1.1.1"

    // help message
    if (params.help) {
        printHelpMessage(NEMU_VERSION, params)
        System.exit(0)
    }

    // Dependency Checks
    def reqs = ["seqkit", "taxonkit", "tblastn", "blastdbcmd", "mafft", "macse",
                "goalign", "python3", "java", "run_treeshrink.py", 
                "nw_reroot", "nw_distance", "nw_prune", "iqtree2", "jq",
                "collect_mutations.py", "calculate_mutspec.py"]
    
    reqs.each { dep ->
        if (!commandExists(dep)) {
            log.error "Required dependency '${dep}' not found in PATH."
            System.exit(1)
        }
    }

    if (params.inputType == "protein" || params.inputType == "prot") {

        log.info """\
            N E M U   P I P E L I N E  ${NEMU_VERSION}
            ================================
            input type   : ${params.inputType}
            input file   : ${params.input}
            blast db     : ${params.db}
            taxdump      : ${params.taxdump}
            max targets  : ${params.maxTargetSeqs}
            gencode      : ${params.gencode}
            MSA mode     : ${params.msaMode}
            IQ-TREE model: ${params.model}
            ASR model    : ${params.modelAsr}
            Threads      : ${params.threads}
            """
            .stripIndent()

        def combined = ["Input file": params.input,
                        "BLAST database": params.db + ".ndb", // TODO must be absolute path
                        "Taxdump directory": params.taxdump]
        combined.each { label, path ->
            if (!path || path == "" || !file(path).exists()) {
                log.error "${label} path '${path}' is empty or does not exist."
                System.exit(1)
            }
        }
        if (!params.msaMode || !(params.msaMode in ["auto", "macse", "mafft_macse", "mafft"])) {
            log.error "Invalid MSA mode specified. Set --msa-mode to 'auto', 'macse', 'mafft_macse', or 'mafft'."
            System.exit(1)
        }

        // Blast + Filter + Extract Nucleotide Sequences
        blastHead(
            params.input, params.speciesName, 
            params.taxdump, params.db, params.maxTargetSeqs, 
            params.gencode
        )

        nuc_multifasta = blastHead.out.nuc_multifasta
        parsed_taxonomy = blastHead.out.taxonomy
        treefile = ""

    } else if (params.inputType == "nucleotide_coding" || 
               params.inputType == "nucleotide_noncoding" ||
               params.inputType == "cds" || params.inputType == "CDS" ||
               params.inputType == "noncds" || params.inputType == "NONCDS"
        ) {
        log.info """\
        N E M U   P I P E L I N E  ${NEMU_VERSION}
        ================================
        input type   : ${params.inputType}
        input file   : ${params.input}
        gencode      : ${params.gencode}
        aligned      : ${params.aligned}
        MSA mode     : ${params.msaMode}
        IQ-TREE model: ${params.model}
        ASR model    : ${params.modelAsr}
        Threads      : ${params.threads}
        """
        .stripIndent()

        treefile = params.treefile
        parsed_taxonomy = null

        input_fasta = channel.fromPath(params.input) // TODO check existance of input files (currently there is no check)
            .filter { fasta -> 
            if (fasta.countFasta() > params.minSeqs) return true
            log.warn "Input file ${fasta.getName()} has less than ${params.minSeqs} sequences. SKIPPING."
            return false
        }
        CHECK_INPUT_TYPE(input_fasta)
        input_fasta_nuc = CHECK_INPUT_TYPE.out.filter { 
            fasta, type ->
            if (type == "DNA") return true
            log.warn "Input file ${fasta.getName()} does not appear to be DNA sequences. SKIPPING."
            return false
        }.map { fasta, _type -> fasta }

        nuc_multifasta = input_fasta_nuc.map { file ->
            def name = file.getBaseName().replaceAll(/[^a-zA-Z0-9]/, '_')
            def outgroupId = params.outgroupId

            // Check if the file contains the outgroupId
            def contains_outgroup = file.text.contains(params.outgroupId)
            if (!contains_outgroup) {
                outgroupId = ""
                log.warn "Outgroup ID '${params.outgroupId}' not found in the file ${name}. This may lead to incorrect rooting of the tree and inaccurate mutation spectra. Continue anyway."
            }
            [name, file, outgroupId]
        }
    }
    else {
        log.error "Invalid input type specified. Set --input-type to 'protein', 'cds', or 'noncds'."
        System.exit(1)
    }

    // in case of protein input, alignment is always needed 
    def aligned = (params.inputType == "protein" || params.inputType == "prot") ? false : params.aligned

    // For noncoding nucleotide input, restrict to mafft (non-codon-aware) alignment
    def effective_msa_mode = (params.inputType == "nucleotide_noncoding" || params.inputType == "noncds" || params.inputType == "NONCDS") ? "mafft" : params.msaMode
    if (effective_msa_mode != params.msaMode) {
        log.warn "Input type is 'noncds' (nucleotide_noncoding); overriding --msa-mode '${params.msaMode}' to 'mafft' (codon-aware modes are not supported for noncoding sequences)."
    }

    // NEMU Core Workflow: Alignment, Phylogeny, ASR, Mutation Extraction, Spectra Derivation
    nemuCore(nuc_multifasta, params.gencode, 
             params.minSeqs, aligned, effective_msa_mode, 
             params.model, params.modelAsr, 
             treefile, params.runTreeShrink,
             params.probaArg, params.uncertaintyCoef, 
             params.consCatCutoff,
             params.plot, params.internal, params.terminal,
             params.branchSpectra)

    // generate readme
    WRITE_README()

    // final outputs aggregation
    spectra_total = nemuCore.out.syn_spectrum.collectFile(
        name: 'spectra_total.tsv', keepHeader: true, skip: 1)

    publish:
    parsed_taxonomy = parsed_taxonomy
    spectra_total = spectra_total
    spectra_data = nemuCore.out.spectra_data
    spectra_plots = nemuCore.out.spectra_plots
    mutation_logs = nemuCore.out.mutation_logs
    mutation_data = nemuCore.out.mutation_data
    encoded_headers = nemuCore.out.encoded_headers
    msa_tree = nemuCore.out.msa_tree
    tree_images = nemuCore.out.tree_images
    readme = WRITE_README.out
}

output {
    spectra_total { mode 'copy' }
    parsed_taxonomy { mode 'copy' }
    readme { mode 'copy' }
    spectra_plots {
        path { sample -> "${sample[0]}/images/" }
    }
    tree_images {
        path { sample -> "${sample[0]}/images/" }
    }
    spectra_data {
        path { sample -> "${sample[0]}/" }
    }
    mutation_logs {
        path { sample -> "${sample[0]}/" }
    }
    mutation_data {
        path { sample -> "${sample[0]}/" }
    }
    encoded_headers {
        path { sample -> "${sample[0]}/" }
    }
    msa_tree {
        path { sample -> "${sample[0]}/" }
    }
}
