#!/usr/bin/env nextflow

NEMU_VERSION="1.1.0"

/*
 * Pipeline: Protein Sequence Analysis
 * Steps:
 * 1. Split Multi-FASTA & Assign Unique IDs
 * 2. Validate Protein sequences
 * 3. Parse Species Name
 * 4. Prepare Taxonomy (TaxonKit: Species & Relatives)
 * 5. BLAST Species & Outgroup (Two-step search)
 */

/* Requirements:
 *
 * Nextflow
 * seqkit v2.9.0
 * taxonkit v0.20.0
 * BLAST+ 2.17.0
 * Python 3
 * mafft
 * trimAl
 * iqtree 2.2.0
 */


// Global parameters
params.input        = "data/proteins.fa"
params.outdir       = "results"
params.db           = "/home/kpotoh/.nuc_db/dolphin"
params.taxdump      = "/home/kpotoh/.taxonkit"
// params.scipts       = "/home/kpotoh/nemu-pipeline/scripts/perl" TODO remove scripts and perl dependencies
params.species_name = false  // Override species name if provided (optional)
params.gencode      = 2      
params.max_target_seqs = 2000 // Limit for species hits
params.min_seqs     = 4    // Minimum number of sequences after filtering to proceed to MSA
params.msa_mode     = "accurate" // "accurate" or "fast"
params.threads      = 4

// Internal parameters
thr_gaps = 0.05 // Gap threshold for MSA cleaning


// check params
if (!params.input || params.input == "") {
    log.error "ERROR: Input file not specified. Use --input to provide input FASTA file."
    System.exit(1)
}
if (!params.db || params.db == "") {
    log.error "ERROR: BLAST database path not specified. Use --db to provide BLAST database."
    System.exit(1)
}
if (!params.msa_mode || !(params.msa_mode in ["auto_cdn", "accurate_cdn", "fast_cdn", "pure_mafft"])) {
    log.error "ERROR: Invalid MSA mode specified. Use --msa_mode with 'auto_cdn', 'accurate_cdn', 'fast_cdn', or 'pure_mafft'."
    System.exit(1)
}

log.info """\
    N E M U   P I P E L I N E  ${NEMU_VERSION}
    =================================
    input file   : ${params.input}
    outdir       : ${params.outdir}
    blast db     : ${params.db}
    min seqs     : ${params.min_seqs}
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
    
    # ---------------------------------------------------------
    # 1. Species BLAST
    # ---------------------------------------------------------
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

    # ---------------------------------------------------------
    # 2. Outgroup BLAST (Relatives)
    # ---------------------------------------------------------
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

    # ---------------------------------------------------------
    # 3. Filter & Combine (Python)
    # ---------------------------------------------------------
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
                
                # saccver pident length qlen gapopen sstart send evalue bitscore sframe
                sacc = parts[0]
                pident = float(parts[1])
                length = float(parts[2])
                qlen = float(parts[3])
                bitscore = float(parts[8])
                sframe = int(parts[9])
                sstart = parts[5]
                send = parts[6]
                
                coverage = length / qlen
                
                if pident >= min_ident and coverage >= min_cov:
                    hits_list.append({
                        'sacc': sacc, 'sstart': sstart, 'send': send, 
                        'sframe': sframe, 'bitscore': bitscore, 'pident': pident
                    })
    except FileNotFoundError:
        pass

# Parse Species and Outroup
parse_blast('blast_species.tsv',  hits,  MIN_IDENT_SPECIES,  MIN_COV_SPECIES)
parse_blast('blast_outgroup.tsv', outgroup_hits, MIN_IDENT_OUTGROUP, MIN_COV_OUTGROUP)

if outgroup_hits:
    # Sort Outgroup by bitscore (descending) to find the Closest Relative
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

    # ---------------------------------------------------------
    # 4. Extract sequences
    # ---------------------------------------------------------
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
    # Encode IDs with seqkit https://bioinf.shenwei.me/seqkit/usage/#replace (Rename with number of record)
    seqkit replace -p .+ -r "seq_{nr}" -w 0 < $sequences > encoded_raw.fasta
    if grep -q "Selected Outgroup" $flt_log; then
        echo "Outgroup sequence included."
        outgrp_id=\$(seqkit seq -i -n < ./encoded_raw.fasta | tail -1)
        seqkit replace -p \${outgrp_id} -r "OUTGRP" -w 0 < encoded_raw.fasta > encoded.fasta
    else
        echo "No outgroup sequence included."
        mv encoded_raw.fasta encoded.fasta
    fi

    # Save mapping of encoded headers TODO
    codes=\$(seqkit seq -ni < ./encoded.fasta)
    original_ids=\$(seqkit seq -ni < ./${sequences})
    paste <(echo "\$codes") <(echo "\$original_names") > encoded_headers.txt
    
    # Remove duplicates
    seqkit rmdup -D duplicated.txt -s -w 0 < encoded.fasta > seqs_unique.fasta
    if grep -q OUTGRP ./duplicated.txt; then
        echo "WARNING: Outgroup sequence was duplicated and removed."
    fi

    # check number of sequences
    num_seqs=\$(grep -c '>' ./seqs_unique.fasta)
    """
}

process MSA_DUMMY {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'

    input:
    tuple val(id), path(sequences), path(encoded_headers), val(num_seqs)
    val gencode
    val msa_mode
    val thr_gaps

    output:
    tuple val(id), path("msa_nuc.fasta")

    script:
    """
    echo "Dummy MSA process for ${id}" > msa_nuc.fasta
    """
}

process MSA {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'
    cpus params.threads

    input:
    tuple val(id), path(sequences), path(encoded_headers), val(num_seqs)
    val gencode
    val msa_mode
    val thr_gaps

    output:
    tuple val(id), path("msa_nuc.fasta"), path("")
    //  file "seq_dd_AA.fa" into seq_dd_AA_for_QC
    //  file "*.csv" optional true

    script:
    """
    if [ $msa_mode = "auto_cdn" ]; then
        if [ ${num_seqs} -le 100 ]; then
            msa_mode_sh="accurate_cdn"
        else
            msa_mode_sh="fast_cdn"
        fi
    else
        msa_mode_sh="$msa_mode"
    fi

    if [[ \$msa_mode_sh = "pure_mafft" ]]; then
        mafft --thread ${task.cpus} $sequences > msa_raw.fasta

    elif [[ \$msa_mode_sh = "accurate_cdn" ]]; then
        mafft --thread ${task.cpus} $sequences > seqM.fa
        sed '/^>/!s/[actg]/\\U&/g' seqM.fa > seqMU.fa
        
        # TODO move filters to the end of the process

        goalign clean seqs -c 0.3 -i seqMU.fa -o seqMC.fa
        goalign clean sites -c $thr_gaps -i seqMC.fa -o seqMCC.fa
        seqkit rmdup -s < seqMCC.fa > seq_dd.fa
        java -jar /opt/macse_v2.07.jar -prog alignSequences -seq seq_dd.fa \
            -gc_def $gencode -optim 2 -max_refine_iter 0 -ambi_OFF
        
        java -jar /opt/macse_v2.07.jar -prog exportAlignment \
            -align seq_dd_NT.fa -gc_def $gencode -ambi_OFF \
            -codonForInternalStop "NNN" -codonForFinalStop "---" \
            -codonForInternalFS "NNN" -codonForExternalFS "---" \
            -out_stat_per_seq macse_stat_per_seq.csv -out_stat_per_site macse_stat_per_site.csv 

        sed 's/!/n/g' seq_dd_NT_NT.fa > seq_dd_NT_FS.fa
        goalign clean sites -c $thr_gaps -i seq_dd_NT_FS.fa -o seq_dd_NT_FS_clean.fa
        seqkit rmdup -s < seq_dd_NT_FS_clean.fa > msa_nuc_lower.fasta
        sed '/^>/!s/[actg]/\\U&/g' msa_nuc_lower.fasta > msa_nuc.fasta

    elif [[ \$msa_mode_sh = "fast_legacy" ]]; then

        # TODO improve according to recomentations of macse team https://www.agap-ge2pop.org/reportgapsaa2nt/

        # NT2AA
        java -jar /opt/macse_v2.07.jar -prog translateNT2AA -seq $sequences \
            -gc_def $gencode -out_AA translated.faa
        #ALN AA
        mafft --thread ${task.cpus} translated.faa > translated_aln.faa
        #AA_ALN --> NT_ALN
        java -jar /opt/macse_v2.07.jar -prog reportGapsAA2NT \
            -align_AA translated_aln.faa -seq $sequences -out_NT aln.fasta
        echo "Do quality control" >&2
        /opt/scripts_latest/macse2.pl aln.fasta msa_nuc.fasta

        cp translated_aln.faa seq_dd_AA.fa
    elif [[ \$msa_mode_sh = "fast_cdn" ]]; then
        # "Big Data" workflow

        MACSE_JAR="/opt/macse_v2.07.jar"

        # TODO add gencode support

        # 1. Trim junk (Optional but recommended)
        java -jar $MACSE_JAR -prog trimNonHomologousFragments \
            -seq $sequences -out_NT 1_trimmed.fasta

        # 2. Repair Frameshifts (The "Fast" MACSE run)
        # We use optimization 0 to make it fast; we just want the FS detection.
        java -jar $MACSE_JAR -prog alignSequences \
            -seq 1_trimmed.fasta -out_NT 2_repaired_with_gaps.fasta \
            -max_refine_iter 0 -local_realign_init 0

        # 3. Remove Gaps (Preserve frameshift corrections)
        # We remove dashes (-) but keep the sequence.
        sed 's/-//g' 2_repaired_with_gaps.fasta > 3_repaired_ungapped.fasta

        # 4. Translate (Now safe because FS are fixed)
        java -jar $MACSE_JAR -prog translateNT2AA \
            -seq 3_repaired_ungapped.fasta -out_AA 4_translated.faa

        # 5. Fast Alignment (MAFFT)
        mafft --auto --thread 4 4_translated.faa > 5_aligned.faa

        # 6. Back-Translate
        java -jar $MACSE_JAR -prog reportGapsAA2NT \
            -align_AA 5_aligned.faa -seq 3_repaired_ungapped.fasta \
            -out_NT final_codon_alignment.fasta

        # Cleanup intermediate files
        rm 1_trimmed.fasta 2_repaired_with_gaps.fasta 3_repaired_ungapped.fasta 4_translated.faa
    
    fi

    # Final cleaning TODO
    #goalign clean seqs -c 0.3 -i seqMU.fa -o seqMC.fa
    #goalign clean sites -c $thr_gaps -i seqMC.fa -o seqMCC.fa
    #seqkit rmdup -s < seqMCC.fa > seq_dd.fa

    """
}


workflow {
    // TODO make 2 workflows: main protein and main nucleotide

    def seq_counter = 0

    raw_sequences = Channel.fromPath(params.input)
        .splitFasta(record: [id: true, header: true, seqString: true])
        .filter { record ->
            def seq = record.seqString.toUpperCase()
            if (seq =~ /[EFILPQZ]/) return true
            def dna_count = seq.count('A') + seq.count('C') + seq.count('G') + seq.count('T') + seq.count('N')
            if (seq.length() > 0 && (dna_count / seq.length()) > 0.95) {
                log.warn "SKIPPING ${record.id}: Sequence appears to be mostly DNA."
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

    // FILTER: Terminate if TaxID is missing
    tax_verified_ch = PREPARE_TAXONOMY.out.filter { id, sp_tax, rel_tax, seq, sp_name ->
        if (sp_tax.size() > 0) {
            return true
        } else {
            log.warn "SKIPPING ${id}: No valid TaxID found for species '${sp_name}'."
            return false
        }
    }

    SAVE_QUERY(tax_verified_ch)
    TBLASTN_AND_FILTER(SAVE_QUERY.out, params.db)

    // FILTER: Terminate if sequences are not found
    seq_verified_ch = TBLASTN_AND_FILTER.out.filter { id, seq, flt_log ->
        if (seq.size() > 0) {
            return true
        } else {
            log.warn "SKIPPING ${id}: No sequences found with TBLASTN for '${id}'."
            return false
        }
    }

    ENCODE_AND_RMDUP(seq_verified_ch)

    // FILTER: Terminate if number of sequences is less than CUTOFF (params.min_seqs)
    seq_num_verified_ch = ENCODE_AND_RMDUP.out.filter { id, seq, enc_head, num_seqs ->
        if (num_seqs.toInteger() > params.min_seqs) {
            return true
        } else {
            log.warn "SKIPPING ${id}: Too low number of sequences for '${id}': ${num_seqs}."
            return false
        }
    }

    MSA_DUMMY(seq_num_verified_ch, params.gencode, params.msa_mode, thr_gaps)
}