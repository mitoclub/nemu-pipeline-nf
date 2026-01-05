#!/usr/bin/env nextflow

/*
 * NeMu: a comprehensive pipeline for accurate reconstruction of neutral mutation spectra from evolutionary data
 * https://doi.org/10.1093/nar/gkae438
 * Authors: Bogdan Efimenko, Konstantin Popadin, Konstantin Gunbin
 */

/* Requirements:
 * Nextflow, seqkit, taxonkit, BLAST+, Python 3.12, mafft, macse, goalign, iqtree2, newick_utils, PyMutSpec
 * Pymutspec (Python lib)
 */

// --- Global Parameters ---

// Inputs/Outputs
params.input            = ""
params.outdir           = "results"

// Databases & Tools
params.db               = ""
params.taxdump          = ""

// Pipeline Logic
params.input_type       = "protein"                     // protein, nucleotide_coding, nucleotide_noncoding
params.species_name     = ""                            // Override species name
params.gencode          = 1
params.max_target_seqs  = 2000
params.min_seqs         = 4                             // Min sequences to proceed
params.threads          = 1

// MSA & Tree
params.outgroup_id      = "OUTGRP"                      // Manually specify outgroup sequence ID for 'nucleotide' input
params.aligned          = false                         // 'nucleotide' input can be pre-aligned
params.msa_mode         = "auto"                        // auto, macse (accurate codon alignment), mafft_macse (fast codon alignment), mafft
params.treefile        = ""                            // Provide pre-computed tree path to skip tree building
params.model            = "GTR+FO+G4+I"                 // IQ-TREE Model
params.model_asr        = "GTR+FO+G4+I"                 // ASR Model
params.run_treeshrink   = true                          // Run TreeShrink to prune long branches

// Mutation Extraction
params.cons_cat_cutoff  = 0                             // 0 = no cutoff // TODO pass list of categories instead of single value
params.proba_arg        = true                          // Use probabilities
params.uncertainty_coef = true                          // Use phylogeny uncertainty coefficient TODO improve implementation

// Spectra Calculation
params.plot             = false                         // Generate plots

// Subsets of mutations to derive spectra for
params.internal         = false
params.terminal         = false
params.branch_spectra   = false


process PARSE_SPECIES_NAME {
    tag "$id"

    input:
    tuple val(id), val(header), val(sequence)
    val global_species

    output:
    tuple val(id), stdout, val(sequence)

    script:
    """
    #!/usr/bin/env python3
    import re
    import sys

    g_spec = "${global_species}"
    if g_spec and g_spec != "false":
        print(g_spec, end='')
        sys.exit(0)

    header = "${header}"
    species = "unknown_species"

    # Try Uniprot style: OS=Homo sapiens
    m_os = re.search(r'OS=([a-zA-Z0-9_ ]+)', header)
    if m_os:
        species = m_os.group(1).strip()
    else:
        # Try NCBI style: [Homo sapiens]
        # Double backslashes are needed for Groovy string interpolation
        m_br = re.search(r'\\[([a-zA-Z0-9_ ]+)\\]', header)
        if m_br:
            species = m_br.group(1).strip()
    
    print(species, end='')
    """
}

process PREPARE_TAXONOMY {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy' // TODO remove line after debug completion

    input:
    tuple val(id), val(species_name), val(sequence)
    path taxdump_dir

    output:
    tuple val(id), path("query.fa"), path("species.taxid"), path("relatives.taxid")

    script:
    """
    if [ "${species_name}" = "unknown_species" ]; then
        touch species.taxid relatives.taxid query.fa
        exit 0
    fi

    export TAXONKIT_DB=${taxdump_dir}
    CLEAN_NAME=\$(echo "${species_name}" | tr '_' ' ')

    echo "Deriving TaxIDs for: \$CLEAN_NAME"
    
    # 1. Get TaxID
    echo "\$CLEAN_NAME" | taxonkit name2taxid | taxonkit lineage -i 2 -r -L > taxid_lineage.txt
    SPEC_ID_RAW=\$(cut -f2 taxid_lineage.txt)
    SPEC_RANK=\$(cut -f3 taxid_lineage.txt)

    if [ -z "\$SPEC_ID_RAW" ]; then
        echo "WARNING: TaxID not found for \$CLEAN_NAME"
        touch species.taxid relatives.taxid query.fa
        exit 0
    fi
    
    if [ \$SPEC_RANK = "species" ]; then
        SPEC_ID=\$SPEC_ID_RAW
    else
        # Get species from lineage
        SPEC_ID=\$(echo \$SPEC_ID_RAW | taxonkit lineage | taxonkit reformat -t -f "{s}" | cut -f4)
        
        # if spec_id is empty, fallback to original
        if [ -z "\$SPEC_ID" ]; then
            SPEC_ID=\$SPEC_ID_RAW
        fi
    fi
    
    # 2. Get downstream TaxIDs
    taxonkit list --ids \$SPEC_ID --indent "" | head -n -1 > species.taxid

    # 3. Find Family ID
    FAMILY_ID=\$(echo \$SPEC_ID | taxonkit lineage | taxonkit reformat -t -f "{f}" | cut -f4)

    if [ -z "\$FAMILY_ID" ]; then
        echo "WARNING: Family rank not found for \$CLEAN_NAME"
        touch relatives.taxid
    else
        # 4. Get Family members and exclude self
        taxonkit list --ids \$FAMILY_ID --indent "" | head -n -1 > family_all.taxid
        grep -vFf species.taxid family_all.taxid > relatives.taxid
    fi

    # Save query
    echo ">${id} ${species_name}" > query.fa
    echo "${sequence}" >> query.fa
    """
}

process TBLASTN_AND_FILTER {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'
    cpus params.threads

    input:
    tuple val(id), path(query), path(species_taxids), path(relatives_taxids)
    val db_path

    output:
    tuple val(id), path("sampled_sequences.fasta"), env("OUTGRP_ID"), emit: seqs
    path "filtering_log.txt"

    script:
    """
    # Define formatting
    outfmt="6 saccver pident length qlen gapopen sstart send evalue bitscore sframe"
    
    # 1. Species BLAST
    if [ -s ${species_taxids} ]; then
        echo "Running Species BLAST..." > filtering_log.txt
        tblastn -query ${query} -db ${db_path} -db_gencode ${params.gencode} \
            -max_target_seqs ${params.max_target_seqs} \
            -evalue 0.0001 -num_threads ${task.cpus} \
            -taxidlist ${species_taxids} -no_taxid_expansion \
            -outfmt "\$outfmt" \
            -out blast_species.tsv
    else
        echo "No species taxids found. Skipping Species BLAST." >> filtering_log.txt
        touch blast_species.tsv
    fi

    # 2. Outgroup BLAST (Relatives)
    if [ -s ${relatives_taxids} ]; then
        echo "Running Outgroup BLAST..." >> filtering_log.txt
        # Only need top 10 hits to find a good outgroup
        tblastn -query ${query} -db ${db_path} -db_gencode ${params.gencode} \
            -max_target_seqs 10 \
            -evalue 0.0001 -num_threads ${task.cpus} \
            -taxidlist ${relatives_taxids} -no_taxid_expansion \
            -outfmt "\$outfmt" \
            -out blast_outgroup.tsv
    else
        echo "No relative taxids found. Skipping Outgroup BLAST." >> filtering_log.txt
        touch blast_outgroup.tsv
    fi

    # 3. Filter & Combine (Python)
    python3 -c "
import sys

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

parse_blast('blast_species.tsv',  hits,  MIN_IDENT_SPECIES,  MIN_COV_SPECIES)
parse_blast('blast_outgroup.tsv', outgroup_hits, MIN_IDENT_OUTGROUP, MIN_COV_OUTGROUP)

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
    " >> filtering_log.txt
    # END OF PYTHON CODE

    OUTGRP_ID=\$(grep -oP "Selected Outgroup: \\K[^,]+" filtering_log.txt)

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
    publishDir "${params.outdir}/${id}", mode: 'copy'

    input:
    tuple val(id), path(sequences), val(OUTGRP_ID)

    output:
    tuple val(id), path("seqs_unique.fasta"), path("encoded_headers.txt"), env('NUM_SEQS')

    script:
    """
    seqkit replace -p .+ -r "seq_{nr}" -w 0 < $sequences > encoded_raw.fasta
    if [ -z "${OUTGRP_ID}" ]; then
        mv encoded_raw.fasta encoded.fasta
    else
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

    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(id), path(sequences)
    val model
    val run_shrinking
    path treefile  // can be error due to automatic path existence check !!!!!! TODO fix this error

    output:
    tuple val(id), path("msa_filtered.fasta"), path("tree.nwk")

    script:
    QUANTILE=0.1
    """
    if [ "${treefile}" != "" ]; then
        echo "Using provided treefile: ${treefile}"
        cp ${treefile} tree.nwk
    else
        echo "Building tree de novo..."
        iqtree2 -s $sequences -m $model -nt $task.cpus --prefix ml

        nseq=\$(grep -c '>' $sequences)
        if [ $run_shrinking = true ] && [ \$nseq -gt 10 ]; then
            run_treeshrink.py -t ml.treefile -O treeshrink -o . -q $QUANTILE -x OUTGRP
            mv treeshrink.treefile treeshrink.nwk
        else
            mv ml.treefile treeshrink.nwk
        fi

        # Check outgroup "quality"
        nw_distance -m p -s f -n treeshrink.nwk | sort -grk 2 > branches.txt
        
        # Prune bad outgroup if needed (simple heuristic: if OUTGRP is not the furthest leaf)
        head -n 1 branches.txt >> branches.head1.txt
        if grep -q OUTGRP branches.head1.txt; then
            nw_reroot -l treeshrink.nwk OUTGRP > tree.nwk
        else
            nw_prune treeshrink.nwk OUTGRP | nw_reroot - > tree.nwk
        fi
    fi

    # drop sequences not present in the tree
    nw_labels -I tree.nwk  > leaves.txt
    seqkit grep -f leaves.txt -w 0 $sequences > msa_filtered.fasta
    seqkit stats $sequences msa_filtered.fasta # sanity check
    """
}

process ASR {
    tag "$id"
    cpus params.threads

    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(id), path(sequences), path(tree)
    val model

    output:
    tuple val(id), path(sequences), path("final_tree.nwk"), path("iqtree_anc.state"), path("rates.tsv")

    script:
    """
    iqtree2 -te $tree -s $sequences -m $model -asr -nt $task.cpus --prefix asr --rate
    mv asr.rate rates.tsv

    if grep -q OUTGRP asr.treefile; then
        nw_reroot -l asr.treefile OUTGRP | sed 's/;/ROOT;/' > final_tree.nwk
    else
        nw_reroot asr.treefile | sed 's/;/ROOT;/' > final_tree.nwk
    fi

    iqtree_states_add_part.py asr.state iqtree_anc.state
    """
}

process DRAW_TREE {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'

    input:
    tuple val(id), path(tree), path(nodes_mapping)

    output:
    path("tree.svg")
    path("tree.png")

    script:
    """
    ## cut -f2 $nodes_mapping  | cut -d ' ' -f 2-3 | sed 's/[^a-zA-Z0-9\\_]/_/g' | cut -c 1-20 > cleaned_headers.txt
    
    cut -f1 -d: $nodes_mapping | sed 's/\t/_/' | cut -c 1-20 > cleaned_headers.txt
    paste <(cut -f1 $nodes_mapping ) <(cat cleaned_headers.txt) > nodes_mapping_cleaned.txt
    
    nw_rename $tree nodes_mapping_cleaned.txt > tree_renamed.nwk
    nw_display -s -b 'visibility:hidden' -i 'visibility:hidden' tree_renamed.nwk > tree.svg
    
    magick tree.svg tree.png # TODO replace by python code (must be less than 80MB of ImageMagick)
    """
}

process MUT_EXTRACTION {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'
    cpus params.threads

    input:
    tuple val(id), path(sequences), path(tree), path(internal_states), path(rates)
    val gencode
    val proba_arg
    val uncertainty_coef
    val cons_cat_cutoff 

    output:
    tuple val(id), path("observed_mutations.tsv"), path("expected_freqs.tsv"), emit: mutations
    tuple val(id), path(sequences), path(tree), path("mut_extraction.log")
    path "expected_mutations.tsv.gz"

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
    publishDir "${params.outdir}/${id}", mode: 'copy'

    input:
    tuple val(id), path(obs_muts), path(exp_freqs)
    val plot
    val internal
    val terminal
    val branch_spectra

    output:
    path "*.tsv"
    path "*.pdf", optional: true

    script:
    """
    nmuts=`cat $obs_muts | wc -l`
    if [ \$nmuts -lt 2 ]; then
        echo "ERROR: There are no reconstructed mutations." >&2
        exit 1 # TODO handle this better check line num in the nextflow level
    fi

    ARGS="--exclude OUTGRP,ROOT --mnum192 16 --proba_cutoff 0.3 --syn --syn4f --all --nonsyn"
    if [ $plot = true ]; then
        ARGS="\$ARGS --plot -x pdf"
    fi

    # TODO replace mean_expected_mutations.tsv with exp_muts if needed
    # TODO replace by pure python code
    
    # Main Calculation
    calculate_mutspec.py -b $obs_muts -e $exp_freqs -o . \$ARGS

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
    """
}

// process AMINO_ACID_STUFF TODO

process AGGREGATE_OUTPUTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path 'spectrum'

    output:
    path("spectra_total.tsv")

    script:
    """
    echo "Aggregating final outputs..."
    python3 -c "
import pandas as pd
import glob
all_spectra = []
for file in glob.glob('spectrum*'):
    df = pd.read_csv(file, sep='\\t')
    #species_name = file.split('/')[-2]  # Assuming structure: outdir/species/ms12syn.tsv TODO fix
    species_name = file
    all_spectra.append(df.assign(Species=species_name))
final_df = pd.concat(all_spectra, ignore_index=True)
final_df.to_csv('spectra_total.tsv', sep='\\t', index=False)
    "
    """
}

process CHECK_INPUT_TYPE {
    input:
    path fasta

    output:
    tuple path(fasta), env("TYPE")

    script:
    """
    TYPE=\$(seqkit stats $fasta -T | tail -1 | cut -f3)
    """
}

process WRITE_README {
    publishDir "${params.outdir}", mode: 'copy'

    output:
    path("readme.txt")

    script:
"""
cat > readme.txt <<- EOM
Output structure:

TODO update after all changes

.
├── final_tree.nwk						# Final phylogenetic tree
├── seqs_unique.fasta					# Filtered orthologous sequences
├── msa_nuc.fasta						# Verified multiple sequence alignment
├── headers_mapping.txt					# Encoded headers of sequences
├── encoded_headers.txt					# Encoded headers of sequences (v2 for different versions of input)
├── logs/
│   ├── report.blast					# Tblastn output during orthologs search
│   ├── *.taxids						# Taxids used in taxa-specific blasing in nt; relatives.taxids contains 
│	│									# 	other species from the genus of query and used for outgroup selection
│   ├── iqtree.log						# IQ-TREE logs during phylogenetic tree inference
│   ├── iqtree_report.log				# IQ-TREE report during phylogenetic tree inference
│   ├── iqtree_treeshrink.log			# TreeShrink logs
│   ├── iqtree_pruned_nodes.log			# Nodes pruned from tree by TreeShrink
│   ├── iqtree_anc.log					# IQ-TREE logs during ancestral reconstrution
│   ├── iqtree_anc_report.log			# IQ-TREE report during ancestral reconstrution
│   ├── iqtree_mut_extraction.log		# Logs during mutation extraction process
│   └── branches.txt					# Tree branch lenghts
├── figures
│   ├── ms12syn.pdf						# Barplot with  12-component spectrum on synonymous mutations
│   └── ms192syn.pdf					# Barplot with 192-component spectrum on synonymous mutations
├── tables
│   ├── rates.tsv						# Site rates categories for an alignment
│   ├── expected_freqs.tsv				# Frequencies of substitutions for each tree node genome
│   ├── mean_expexted_mutations.tsv		# Averaged frequencies of substitutions for entire tree
│   ├── ms12syn.tsv						# table with 12-component spectrum on synonymous mutations
│   ├── ms192syn.tsv					# table with 192-component spectrum on synonymous mutations
│   └── observed_mutations.tsv			# Recontructed mutations
EOM
"""
}

// Helper: Check dependencies
boolean commandExists(String command) {
    def proc = ["bash", "-c", "command -v ${command}"].execute()
    proc.waitFor()
    return proc.exitValue() == 0
}

workflow {
       
    NEMU_VERSION="1.1.0"

    // help message
    if (params.help) {
        println """
N E M U   P I P E L I N E  ${NEMU_VERSION}
================================
A comprehensive pipeline for accurate reconstruction of neutral mutation spectra from evolutionary data.
https://doi.org/10.1093/nar/gkae438

Usage: nextflow run main.nf --input <input_fasta> [options]

Main options:
    --input_type STRING     Type of input sequences: 
                            protein, nucleotide_coding, nucleotide_noncoding (default: ${params.input_type})
    --input FILE            Input FASTA file (required)
                            If input type is protein, a multi-FASTA with one or several sequences 
                            required (header format: ">ID [Species name]"). If input type 
                            is nucleotide, one or several fasta files with orthologous 
                            sequences (including outgroup) required 
    --outdir DIR            Output directory (default: ${params.outdir})
    --gencode NUM           Genetic code table (default: ${params.gencode})
                            Used for codon-aware alignment and annotation of mutations

Common options:
    --threads NUM           Number of threads to use (default: ${params.threads})
    --help                  Show this help message and exit

Options for protein input:
    --db PATH               BLAST database path (required for protein input)
    --taxdump DIR           Taxdump directory path (required for protein input)
                            TODO integrate to the container
    --species_name STRING   Override species name. Useful when you work with proteins 
                            from single species
    --max_target_seqs NUM   Max target sequences for BLAST (default: ${params.max_target_seqs})
                            tblastn will collect no more than this number of sequences
    
Options for nucleotide input:
    --outgroup_id STRING    Outgroup sequence ID (default: ${params.outgroup_id})
                            Specify outgroup sequence ID in the alignment for rooting the tree
    --aligned BOOL          Input sequences are pre-aligned (default: ${params.aligned})

Options for MSA & Phylogeny:
    --msa_mode STRING       MSA mode: auto, macse, mafft_macse, mafft (default: ${params.msa_mode})
    --min_seqs NUM          Minimum number of sequences to proceed phylogenetic inference (default: ${params.min_seqs})
    --treefile FILE         Input tree file (optional; default: build tree de novo)
    --model STRING          IQ-TREE substitution model (default: ${params.model})
    --model_asr STRING      ASR substitution model (default: ${params.model_asr})
    --run_treeshrink BOOL   Run TreeShrink to prune long branches (default: ${params.run_treeshrink})

Options for Mutation Extraction:
    --cons_cat_cutoff NUM   Conservation category cutoff for mutation extraction (default: ${params.cons_cat_cutoff})
                            0 = no cutoff; only mutations in sites with rate category 
                            less than or equal to this value will be used
    --proba_arg BOOL        Use probabilities in mutation extraction (default: ${params.proba_arg})
    --uncertainty_coef BOOL Use phylogeny uncertainty coefficient in mutation extraction 
                            (default: ${params.uncertainty_coef})

Options for Mutation Spectra Derivation:
    --plot BOOL             Generate barcharts of mutation spectra (default: ${params.plot})
    --internal BOOL         Derive spectra for internal branches (default: ${params.internal})
    --terminal BOOL         Derive spectra for terminal branches (default: ${params.terminal})
    --branch_spectra BOOL   Derive spectra for individual branches (default: ${params.branch_spectra})
        """.stripIndent()
        System.exit(0)
    }

    // Dependency Checks
    def reqs = ["seqkit", "taxonkit", "tblastn", "blastdbcmd", "mafft", "macse",
                "goalign", "python3", "java", "run_treeshrink.py", 
                "nw_reroot", "nw_distance", "nw_prune", "iqtree2", 
                "collect_mutations.py", "calculate_mutspec.py"]
    
    reqs.each { dep ->
        if (!commandExists(dep)) {
            log.error "Required dependency '${dep}' not found in PATH."
            System.exit(1)
        }
    }

    def combined = ["Input file": params.input,
                    "BLAST database": params.db + ".ndb",
                    "Taxdump directory": params.taxdump]
    combined.each { label, path ->
        if (!path || path == "" || !file(path).exists()) {
            log.error "${label} path '${path}' is empty or does not exist."
            System.exit(1)
        }
    }

    if (params.input_type == "protein") {

    log.info """\
        N E M U   P I P E L I N E  ${NEMU_VERSION}
        ================================
        input type   : ${params.input_type}
        input file   : ${params.input}
        outdir       : ${params.outdir}
        blast db     : ${params.db}
        taxdump      : ${params.taxdump}
        max targets  : ${params.max_target_seqs}
        gencode      : ${params.gencode}
        MSA mode     : ${params.msa_mode}
        IQ-TREE model: ${params.model}
        ASR model    : ${params.model_asr}
        Threads      : ${params.threads}
        """
        .stripIndent()

    if (!params.db || params.db == "") {
        log.error "BLAST database path not specified. Use --db to provide BLAST database."
        System.exit(1)
    }
    if (!params.taxdump || params.taxdump == "") {
        log.error "Taxdump directory not specified. Use --taxdump to provide taxdump path."
        System.exit(1)
    }
    if (!params.msa_mode || !(params.msa_mode in ["auto", "macse", "mafft_macse", "mafft"])) {
        log.error "Invalid MSA mode specified. Use --msa_mode with 'auto', 'macse', 'mafft_macse', or 'mafft'."
        System.exit(1)
    }

    def seq_counter = 0
    // TODO add species name to id ???
    raw_sequences = Channel.fromPath(params.input)
        .splitFasta(record: [id: true, header: true, seqString: true])
        .filter { record ->
            def seq = record.seqString.toUpperCase()
            if (seq =~ /[EFILPQZ]/) return true
            def dna_count = seq.count('A') + seq.count('C') + seq.count('G') + seq.count('T') + seq.count('N')
            if (seq.length() > 0 && (dna_count / seq.length()) > 0.95) {
                log.warn "SKIPPING ${record.id}: Looks like DNA."
                return false
            }
            return true
        }
        .map { record ->
            def count = ++seq_counter
            def clean_original_id = record.id.split()[0].replaceAll(/[^a-zA-Z0-9\.]/, '_')            
            def unique_id = "${count}__${clean_original_id}"
            [unique_id, record.header, record.seqString]
        }

    PARSE_SPECIES_NAME(raw_sequences, params.species_name)
    PREPARE_TAXONOMY(PARSE_SPECIES_NAME.out, params.taxdump)

    // Filter missing TaxIDs
    tax_verified_ch = PREPARE_TAXONOMY.out.filter { id, _q, sp_tax, _rel_tax ->
        if (sp_tax.size() > 0) return true
        log.warn "SKIPPING ${id}: No valid TaxID found."
        return false
    }

    TBLASTN_AND_FILTER(tax_verified_ch, params.db)

    // Filter missing sequences
    fasta_verified_ch = TBLASTN_AND_FILTER.out.seqs.filter { id, fasta, outgrp_id ->
        if (outgrp_id == '') log.warn "WARNING: Outgroup not found for ${id}"
        if (fasta != null && fasta.size() > 0) return true
        log.warn "SKIPPING ${id}: No sequences found."
        return false
    }

    } else if (params.input_type == "nucleotide_coding" || params.input_type == "nucleotide_noncoding") {
    
        log.info """\
        N E M U   P I P E L I N E  ${NEMU_VERSION}
        ================================
        input type   : ${params.input_type}
        input file   : ${params.input}
        outdir       : ${params.outdir}
        gencode      : ${params.gencode}
        aligned      : ${params.aligned}
        MSA mode     : ${params.msa_mode}
        IQ-TREE model: ${params.model}
        ASR model    : ${params.model_asr}
        Threads      : ${params.threads}
        """
        .stripIndent()

        def seq_counter = 0
        input_fasta = Channel.fromPath(params.input)
        input_fasta_nuc = CHECK_INPUT_TYPE(input_fasta).out.filter { 
            fasta, type ->
            if (type == "DNA") return true
            log.warn "Input file ${fasta.getName()} does not appear to be DNA sequences. SKIPPING."
            return false
        }.map { fasta, _type -> fasta }

        fasta_verified_ch = input_fasta_nuc.map { file ->
            def name = file.getBaseName().replaceAll(/[^a-zA-Z0-9]/, '_')
            def count = ++seq_counter
            def name_indexed = "${count}__${name}"

            // Check if the file contains the outgroup_id
            def contains_outgroup = file.text.contains(params.outgroup_id)
            if (!contains_outgroup) {
                log.warn "Outgroup ID '${params.outgroup_id}' not found in the file ${name}. Continue anyway."
            }
            [name_indexed, file, params.outgroup_id]
        }
    }
    else {
        log.error "Invalid input type specified. Use --input_type with 'protein', 'nucleotide_coding' or 'nucleotide_noncoding'."
        System.exit(1)
    }

    ENCODE_AND_RMDUP(fasta_verified_ch)

    // Filter Low Count
    seq_num_verified_ch = ENCODE_AND_RMDUP.out.filter { id, _seq, _enc_head, num_seqs ->
        if (num_seqs.toInteger() > params.min_seqs) return true
        log.warn "SKIPPING ${id}: Count ${num_seqs} < ${params.min_seqs}"
        return false
    }

    // in case of protein input, alignment is always needed
    def aligned = (params.input_type == "protein") ? false : params.aligned

    if (aligned == true) {
        msa_num_verified_ch = seq_num_verified_ch.map { id, seq, _enc_head, _num_seqs ->
            [id, seq]
        }
    } else {
        MSA(seq_num_verified_ch, params.gencode, params.msa_mode)

        msa_num_verified_ch = MSA.out.filter { id, _seq, _msa_log, num_seqs ->
            if (num_seqs.toInteger() > params.min_seqs) return true
            log.warn "SKIPPING ${id}: Count after MSA ${num_seqs} < ${params.min_seqs}"
            return false
        }.map { id, seq, _msa_log, _num_seqs -> [id, seq] }
    }

    WRITE_README()

    // Build or use user-provided tree
    treefile = ""
    if (params.treefile && params.treefile != "") {
        if (!file(params.treefile).exists()) {
            log.warn "User-provided tree file '${params.treefile}' does not exist. Ignoring and building tree de novo."
        } else {
            log.info "Using user-provided tree file: ${params.treefile}"
            treefile = params.treefile
        }
    }
    BUILD_TREE(msa_num_verified_ch, params.model, params.run_treeshrink, treefile)

    ASR(BUILD_TREE.out, params.model_asr)

    // make channel with tree from ASR and nodes mapping from ENCODE_AND_RMDUP for drawing
    just_tree_ch = ASR.out.map { id, _seq, tree, _anc_state, _rates -> [id, tree] }
                          .join(seq_num_verified_ch.map { id, _seq, enc_head, _num_seqs -> [id, enc_head] })
    DRAW_TREE(just_tree_ch)

    MUT_EXTRACTION(ASR.out, 
        params.gencode, params.proba_arg, params.uncertainty_coef,
        params.cons_cat_cutoff,
    )

    DERIVE_SPECTRA(MUT_EXTRACTION.out.mutations, params.plot, 
        params.internal, params.terminal, params.branch_spectra
    )

    // // final outputs aggregation TODO move to separate workflow??
    // spectra = Channel.fromPath( "${params.outdir}/*/ms12syn.tsv" )
    // AGGREGATE_OUTPUTS(spectra)
}