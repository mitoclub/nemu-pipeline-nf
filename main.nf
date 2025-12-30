#!/usr/bin/env nextflow

NEMU_VERSION="1.1.0"

/*
 * Pipeline: Protein Sequence Analysis (NEMU)
 */

/* Requirements:
 * Nextflow, seqkit, taxonkit, BLAST+, Python 3.12, mafft, goalign, iqtree2, newick_utils
 * Pymutspec (Python lib)
 */

// --- Global Parameters ---

// Inputs/Outputs
params.input            = "data/proteins.fa"
params.outdir           = "results"

// Databases & Tools
params.db               = "${System.getenv('HOME')}/.nuc_db/dolphin"
params.taxdump          = "${System.getenv('HOME')}/.taxonkit"
params.macse_jar        = "/opt/macse_v2.07.jar"        // Path to MACSE jar

// Pipeline Logic
params.species_name     = false                         // Override species name
params.gencode          = 2
params.max_target_seqs  = 2000
params.min_seqs         = 4                             // Min sequences to proceed
params.threads          = 4

// MSA & Tree
params.msa_mode         = "auto_codon"                  // auto_codon, accurate_codon, fast_codon, pure_mafft
params.model            = "GTR+FO+G6+I"                 // IQ-TREE Model
params.model_asr        = "auto"                        // ASR Model (auto = same as tree)
params.run_treeshrink   = true                          // Run TreeShrink to prune long branches

// Mutation Extraction
params.cons_cat_cutoff  = 0                             // 0 = no cutoff // TODO pass list of categories instead of single value
params.proba_arg        = true                          // Use probabilities
params.uncertainty_coef = true                          // Use phylogeny uncertainty coefficient TODO improve implementation
params.save_exp_mutations = false

// Spectra Calculation (DERIVE_SPECTRA)
params.plot             = false                         // Generate plots

// Subsets to calculate
params.internal         = false
params.terminal         = false
params.branch_spectra   = false


if (params.model_asr == "auto") {
    params.model_asr = params.model
}


log.info """\
    N E M U   P I P E L I N E  ${NEMU_VERSION}
    =================================
    input file   : ${params.input}
    outdir       : ${params.outdir}
    blast db     : ${params.db}
    macse jar    : ${params.macse_jar}
    min seqs     : ${params.min_seqs}
    gencode      : ${params.gencode}
    MSA mode     : ${params.msa_mode}
    IQ-TREE model: ${params.model}
    Threads      : ${params.threads}
    """
    .stripIndent()

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
    if g_spec and g_spec != "false" and g_spec != "":
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
    publishDir "${params.outdir}/${id}", mode: 'copy'

    input:
    tuple val(id), val(species_name), val(sequence)
    path taxdump_dir

    output:
    tuple val(id), path("species.taxid"), path("relatives.taxid"), val(sequence), val(species_name)

    script:
    """
    if [ "${species_name}" = "unknown_species" ]; then
        touch species.taxid relatives.taxid
        exit 0
    fi

    export TAXONKIT_DB=${taxdump_dir}
    CLEAN_NAME=\$(echo "${species_name}" | tr '_' ' ')

    echo "Deriving TaxIDs for: \$CLEAN_NAME"
    
    # 1. Get TaxID
    SPEC_ID=\$(echo "\$CLEAN_NAME" | taxonkit name2taxid | cut -f2)
    
    # TODO add here check that SPEC_ID is species taxid if it's subspecies
    # TODO make sure that same species name will have no bugs

    if [ -z "\$SPEC_ID" ]; then
        echo "WARNING: TaxID not found for \$CLEAN_NAME"
        touch species.taxid relatives.taxid
    else
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
    fi
    """
}

process SAVE_QUERY {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'

    input:
    tuple val(id), path(species_taxid), path(relatives_taxid), val(sequence), val(species_name)

    output:
    tuple val(id), path("query.fa"), path("species.txt"), path(species_taxid), path(relatives_taxid)

    script:
    """
    echo ">${id}" > query.fa
    echo "${sequence}" >> query.fa
    echo "${species_name}" > species.txt
    """
}

process TBLASTN_AND_FILTER {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'
    cpus params.threads

    input:
    tuple val(id), path(query), path(species_txt), path(species_taxids), path(relatives_taxids)
    val db_path

    output:
    tuple val(id), path("sampled_sequences.fasta"), path("filtering_log.txt")

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
    tuple val(id), path(sequences), path(flt_log)

    output:
    tuple val(id), path("seqs_unique.fasta"), path("encoded_headers.txt"), env(num_seqs)

    script:
    """
    seqkit replace -p .+ -r "seq_{nr}" -w 0 < $sequences > encoded_raw.fasta
    if grep -q "Selected Outgroup" $flt_log; then
        outgrp_id=\$(seqkit seq -i -n < ./encoded_raw.fasta | tail -1)
        seqkit replace -p \${outgrp_id} -r "OUTGRP" -w 0 < encoded_raw.fasta > encoded.fasta
    else
        mv encoded_raw.fasta encoded.fasta
    fi

    # Save mapping
    codes=\$(seqkit seq -ni < ./encoded.fasta)
    original_names=\$(seqkit seq -n < ./${sequences})
    paste <(echo "\$codes") <(echo "\$original_names") > encoded_headers.txt
    
    # Remove duplicates
    seqkit rmdup -D duplicated.txt -s -w 0 < encoded.fasta > seqs_unique.fasta
    num_seqs=\$(grep -c '>' ./seqs_unique.fasta)
    """
}

// THRESHOLDS
LARGE_DATA_CUTOFF=250    // Switch to Big Data workflow if seqs > this
MIN_SEQ_LEN=100          // Pre-filter: remove sequences shorter than 100bp
MAX_GAP_SEQ=0.50         // Post-filter: Remove SEQS with >50% gaps
MAX_GAP_SITE=0.50        // Post-filter: Remove SITES (columns) with >50% gaps

aln_logfile = "alignment.log"

process MSA {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'
    cpus params.threads

    input:
    tuple val(id), path(sequences), path(encoded_headers), val(num_seqs)
    val gencode
    val msa_mode

    output:
    tuple val(id), path("msa_nuc.fasta"), path(aln_logfile), val(new_num_seqs)

    script:
    """
    echo "--- STARTING MSA ---" > $aln_logfile
    echo "[0/5] Pre-filtering short sequences (<${MIN_SEQ_LEN}bp)..." >> $aln_logfile
    seqkit seq -m $MIN_SEQ_LEN -g "$sequences" > input_clean.fasta
    SEQ_COUNT=\$(grep -c "^>" input_clean.fasta)
    echo "Sequences remaining: \$SEQ_COUNT" >> $aln_logfile

    # MSA Mode Selection
    if [ $msa_mode = "auto_codon" ]; then
        if [ "\$SEQ_COUNT" -gt "$LARGE_DATA_CUTOFF" ]; then
            msa_mode_sh="fast_codon"
        else
            msa_mode_sh="accurate_codon"
        fi
    else
        msa_mode_sh="$msa_mode"
    fi

    # ALIGNMENT
    if [ \$msa_mode_sh = "fast_codon" ]; then
        # Big Data Strategy (Trim -> MacseRepair -> Mafft -> MacseBackTrans)
        java -jar "${params.macse_jar}" -prog trimNonHomologousFragments \
            -seq input_clean.fasta -gc_def "$gencode" \
            -out_NT 1_trimmed.fasta > /dev/null 2>&1

        java -jar "${params.macse_jar}" -prog alignSequences \
            -seq 1_trimmed.fasta -gc_def "$gencode" \
            -out_NT 2_repaired.fasta \
            -max_refine_iter 0 -local_realign_init 0 > /dev/null 2>&1

        sed 's/-//g' 2_repaired.fasta > 3_ungapped.fasta

        java -jar "${params.macse_jar}" -prog translateNT2AA \
            -seq 3_ungapped.fasta -gc_def "$gencode" \
            -out_AA 4_protein.faa > /dev/null 2>&1

        mafft --thread "${task.cpus}" --auto --quiet 4_protein.faa > 5_aligned_protein.faa

        java -jar "${params.macse_jar}" -prog reportGapsAA2NT \
            -align_AA 5_aligned_protein.faa \
            -seq 3_ungapped.fasta -gc_def "$gencode" \
            -out_NT raw_alignment.fasta > /dev/null 2>&1

        # Cleanup intermediate files TODO uncomment
        # rm 1_trimmed.fasta 2_repaired.fasta 3_ungapped.fasta 4_protein.faa 5_aligned_protein.faa

    elif [ \$msa_mode_sh = "accurate_codon" ]; then
        # Pure MACSE
        java -jar "${params.macse_jar}" -prog alignSequences \
            -seq input_clean.fasta -gc_def "$gencode" \
            -out_NT raw_alignment.fasta \
            -out_AA raw_alignment_AA.fasta > /dev/null 2>&1
    
    elif [ \$msa_mode_sh = "pure_mafft" ]; then
        mafft --thread "${task.cpus}" --auto --quiet input_clean.fasta > raw_alignment.fasta
    fi

    # SANITIZING
    if [ \$msa_mode_sh != "pure_mafft" ]; then
        java -jar "${params.macse_jar}" -prog exportAlignment \
            -align raw_alignment.fasta \
            -gc_def "$gencode" \
            -ambi_OFF \
            -codonForInternalStop "NNN" -codonForFinalStop "---" \
            -codonForInternalFS "NNN" -codonForExternalFS "---" \
            -out_NT sanitized_alignment.fasta \
            -out_stat_per_seq macse_stat_per_seq.csv \
            -out_stat_per_site macse_stat_per_site.csv > /dev/null 2>&1
    else
        mv raw_alignment.fasta sanitized_alignment.fasta
    fi

    # POST-ALIGNMENT FILTERING
    # Remove columns where >50% of sequences have a gap. 
    # Remove sequences that are >50% gaps (after site cleaning).
    goalign clean sites -c "$MAX_GAP_SITE" -i sanitized_alignment.fasta -o filtered_sites.fasta
    goalign clean seqs -c "$MAX_GAP_SEQ" -i filtered_sites.fasta -o "filtered_seqs.fasta"
    seqkit rmdup -s < filtered_seqs.fasta > msa_nuc.fasta

    # Stats
    echo "--- Raw Alignment Stats ---" >> $aln_logfile
    seqkit stats raw_alignment.fasta >> $aln_logfile
    echo "--- Filtered Alignment Stats ---" >> $aln_logfile
    seqkit stats "msa_nuc.fasta" >> $aln_logfile
    
    # N-Content Warning
    seqkit fx2tab --name --gc --avg-qual "msa_nuc.fasta" | awk '\$4 > 20 {print \$1 " has high N content"}' >> $aln_logfile
    new_num_seqs=\$(grep -c "^>" msa_nuc.fasta)
    """
}

QUANTILE=0.1

process BUILD_TREE {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'
    cpus params.threads

    input:
    tuple val(id), path(sequences)
    val model
    val run_shrinking

    output:
    tuple val(id), path(sequences), path("tree.nwk"), path("iqtree.report")

    errorStrategy 'retry'
    maxRetries 3

    script:
    """
    iqtree2 -s $sequences -m $model -nt $task.cpus --prefix iqtree

    nseq=\$(grep -c '>' $sequences)
    if [ $run_shrinking = true ] && [ \$nseq -gt 10 ]; then
        run_treeshrink.py -t iqtree.treefile -O treeshrink -o . -q $QUANTILE -x OUTGRP
    else
        cat iqtree.treefile > treeshrink.nwk
    fi

    # Check outgroup "quality"
    nw_distance -m p -s f -n treeshrink.nwk | sort -grk 2 > branches.txt
    head -n 1 branches.txt >> branches.head1.txt
    
    # Prune bad outgroup if needed (simple heuristic: if OUTGRP is not the furthest leaf)
    if grep -q OUTGRP branches.head1.txt; then
        nw_reroot -l treeshrink.nwk OUTGRP > tree.nwk
    else
        nw_prune treeshrink.nwk OUTGRP | nw_reroot - > tree.nwk
    fi
    """
}

process ASR {
    tag "$id"
    cpus params.threads

    input:
    tuple val(id), path(sequences), path(tree), path(iqtree_report)
    val model

    output:
    tuple val(id), path('msa_filtered.fasta'), path("final_tree.nwk"), path("iqtree_anc.state"), path("rates.tsv")

    errorStrategy 'retry'
    maxRetries 3

    script:
    """
    # drop sequences not present in the tree
    nw_labels -I $tree | sed 's/\$/\$/' > leaves.txt
    seqkit grep -n -f leaves.txt -w 0 $sequences > msa_filtered.fasta

    iqtree2 -te $tree -s msa_filtered.fasta -m $model -asr -nt $task.cpus --prefix asr --rate
    mv asr.rate rates.tsv

    if grep -q OUTGRP anc.treefile; then
        nw_reroot -l anc.treefile OUTGRP | sed 's/;/ROOT;/' > final_tree.nwk
    else
        nw_reroot anc.treefile | sed 's/;/ROOT;/' > final_tree.nwk
    fi

    iqtree_states_add_part.py anc.state iqtree_anc.state
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
    val save_exp_mutations

    output:
    tuple val(id), path("observed_mutations.tsv"), path("expected_freqs.tsv"), path("mut_extraction.log"), emit: mutations
    tuple val(id), path(sequences), path(tree)
    path "expected_mutations.tsv", optional: true

    script:
    """
    ARGS="--gencode $gencode --no-mutspec --outdir mout --threads ${task.cpus} --syn --syn4f --all --nonsyn"
    if [ $cons_cat_cutoff -gt 0 ]; then
        ARGS="\$ARGS --rates $rates --cat-cutoff $cons_cat_cutoff"
    fi
    if [ $proba_arg = "true" ]; then
        ARGS="\$ARGS --proba --pcutoff 0.3"
    fi
    if [ $save_exp_mutations = "true" ]; then
        ARGS="\$ARGS --save-exp-muts"
    fi
    if [ $uncertainty_coef = "true" ]; then
        ARGS="\$ARGS --phylocoef"
    else
        ARGS="\$ARGS --no-phylocoef"
    fi

    collect_mutations.py --tree $tree --states $sequences --states $internal_states \
        \$ARGS 

    mv mout/* .
    mv mutations.tsv observed_mutations.tsv
    mv run.log mut_extraction.log
    """
}

process DERIVE_SPECTRA {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'

    input:
    tuple val(id), path(obs_muts), path(exp_freqs)
    val plot

    output:
    path "*.tsv"
    path "*.pdf", optional: true

    script:
    """
    nmuts=`cat $obs_muts | wc -l`
    if [ \$nmuts -lt 2 ]; then
        echo "ERROR: There are no reconstructed mutations." >&2
        exit 1
    fi

    # Main Calculation
    calculate_mutspec.py -b $obs_muts -e $exp_freqs -o . \
        --exclude OUTGRP,ROOT --mnum192 16 \
        --proba_cutoff 0.3 --plot -x pdf \
        --syn --syn4f --all --nonsyn

    # Internal
    if [ "$params.internal" = "true" ]; then
        calculate_mutspec.py -b $obs_muts -e $exp_freqs -o . \
            --exclude OUTGRP,ROOT --mnum192 16 \
            --proba_cutoff 0.3 --plot -x pdf --subset internal \
            --syn --syn4f --all --nonsyn
        # Cleanup
        if [ -f mean_expexted_mutations_internal.tsv ]; then rm mean_expexted_mutations_internal.tsv; fi
    fi
    
    # Terminal
    if [ "$params.terminal" = "true" ]; then
        calculate_mutspec.py -b $obs_muts -e $exp_freqs -o . \
            --exclude OUTGRP,ROOT --mnum192 16 \
            --proba_cutoff 0.3 --plot -x pdf --subset terminal \
            --syn --syn4f --all --nonsyn
        if [ -f mean_expexted_mutations_terminal.tsv ]; then rm mean_expexted_mutations_terminal.tsv; fi
    fi
    
    # Branch Spectra
    if [ "$params.branch_spectra" = "true" ]; then
        calculate_mutspec.py -b $obs_muts -e $exp_freqs -o . \
            --exclude OUTGRP,ROOT --mnum192 16 \
            --proba_cutoff 0.3 --branches \
            --syn --syn4f --all --nonsyn
    fi
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
    // TODO make 2 workflows: main protein and main nucleotide
    
    // Dependency Checks
    def reqs = ["seqkit", "taxonkit", "tblastn", "blastdbcmd", "mafft", 
                "goalign", "python3", "java", "run_treeshrink.py", 
                "nw_reroot", "nw_distance", "nw_prune", "iqtree2", 
                "collect_mutations.py", "calculate_mutspec.py"]
    
    for (dep in reqs) {
        if (!commandExists(dep)) {
            log.error "ERROR: Required dependency '${dep}' not found in PATH."
            System.exit(1)
        }
    }

    // Check MACSE JAR file
    if (!file(params.macse_jar).exists()) {
        log.error "ERROR: MACSE JAR file not found at path: ${params.macse_jar}"
        System.exit(1)
    }

    // check params TODO check file existence
    if (!params.input || params.input == "") {
        log.error "ERROR: Input file not specified. Use --input to provide input FASTA file."
        System.exit(1)
    }
    if (!params.db || params.db == "") {
        log.error "ERROR: BLAST database path not specified. Use --db to provide BLAST database."
        System.exit(1)
    }
    if (!params.msa_mode || !(params.msa_mode in ["auto_codon", "accurate_codon", "fast_codon", "pure_mafft"])) {
        log.error "ERROR: Invalid MSA mode specified. Use --msa_mode with 'auto_codon', 'accurate_codon', 'fast_codon', or 'pure_mafft'."
        System.exit(1)
    }

    WRITE_README()
    def seq_counter = 0
    
    // 1. INPUT
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

    // 2. PARSE & TAXONOMY
    PARSE_SPECIES_NAME(raw_sequences, params.species_name)
    PREPARE_TAXONOMY(PARSE_SPECIES_NAME.out, params.taxdump)

    // Filter missing TaxIDs
    tax_verified_ch = PREPARE_TAXONOMY.out.filter { id, sp_tax, rel_tax, seq, sp_name ->
        if (sp_tax.size() > 0) return true
        log.warn "SKIPPING ${id}: No valid TaxID found for species '${sp_name}'."
        return false
    }

    SAVE_QUERY(tax_verified_ch)

    // 3. BLAST
    TBLASTN_AND_FILTER(SAVE_QUERY.out, params.db)

    // Filter missing sequences
    seq_verified_ch = TBLASTN_AND_FILTER.out.filter { id, seq, flt_log ->
        if (seq.size() > 0) return true
        log.warn "SKIPPING ${id}: No sequences found."
        return false
    }

    // 4. CLEANUP
    ENCODE_AND_RMDUP(seq_verified_ch)

    // Filter Low Count
    seq_num_verified_ch = ENCODE_AND_RMDUP.out.filter { id, seq, enc_head, num_seqs ->
        if (num_seqs.toInteger() > params.min_seqs) return true
        log.warn "SKIPPING ${id}: Count ${num_seqs} < ${params.min_seqs}"
        return false
    }

    // 5. MSA
    MSA(seq_num_verified_ch, params.gencode, params.msa_mode)

    msa_num_verified_ch = MSA.out.filter { id, seq, msa_log, num_seqs ->
        if (num_seqs.toInteger() > params.min_seqs) return true
        log.warn "SKIPPING ${id}: Count after MSA ${num_seqs} < ${params.min_seqs}"
        return false
    }.map { id, seq, msa_log, num_seqs -> [id, seq] }

    // 6. TREE & ASR
    BUILD_TREE(msa_num_verified_ch, params.model, params.run_treeshrink)
    ASR(BUILD_TREE.out, params.model_asr)

    // 7. MUTATIONS
    MUT_EXTRACTION(ASR.out, 
        params.gencode, params.proba_arg, params.uncertainty_coef,
        params.cons_cat_cutoff, params.save_exp_mutations, 
    )

    DERIVE_SPECTRA(MUT_EXTRACTION.out.mutations, params.plot)
}