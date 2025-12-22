#!/usr/bin/env nextflow

/*
 * Pipeline: Protein Sequence Analysis
 * Steps:
 * 1. Split Multi-FASTA & Assign Unique IDs
 * 2. Validate Protein sequences
 * 3. Parse Species Name
 * 4. Prepare Taxonomy (TaxonKit: Species & Relatives)
 * 5. BLAST Species & Outgroup (Two-step search)
 */

// Global parameters
params.input        = "data/proteins.fa"
params.outdir       = "results"
params.db           = "/home/kpotoh/.nuc_db/dolphin"
params.taxdump      = "/home/kpotoh/.taxonkit" 
params.species_name = false  
params.gencode      = 2      
params.max_target_seqs = 2000 // Limit for species hits

log.info """\
    P R O T E I N   P I P E L I N E
    ===================================
    input file   : ${params.input}
    outdir       : ${params.outdir}
    blast db     : ${params.db}
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

/*
 * PROCESS: TBLASTN & FILTER
 * 1. Blast Species (limit 2000)
 * 2. Blast Relatives (limit 10)
 * 3. Filter & Extract
 */
process TBLASTN_AND_FILTER {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'
    cpus 4 

    input:
    tuple val(id), path(query), path(species_txt), path(species_taxids), path(relatives_taxids)
    val db_path

    output:
    path "sampled_sequences.fasta"
    path "filtering_log.txt"

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
log_lines = []

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
    log_lines.append(f\\"Selected Outgroup: {best_out['sacc']}, Ident: {best_out['pident']}%, Score: {best_out['bitscore']}\\")
    hits.append(best_out)    
else:
    log_lines.append('WARNING: No valid outgroup found.')

final_coords = []
log_lines = []

log_lines.append(f'Found {len(hits)} valid species hits.')
for h in hits:
    # Logic: ID start-end strand
    strand = 'minus' if h['sframe'] < 0 else 'plus'
    x, y = (h['sstart'], h['send']) if strand == 'plus' else (h['send'], h['sstart'])
    final_coords.append(f\\"{h['sacc']} {x}-{y} {strand}\\")

with open('extract_coords.txt', 'w') as f:
    for line in final_coords:
        f.write(line + '\\n')

with open('filtering_log.txt', 'a') as f:
    f.write('\\n'.join(log_lines) + '\\n')
    "
    # END OF PYTHON CODE

    # ---------------------------------------------------------
    # 4. Extract
    # ---------------------------------------------------------
    if [ -s extract_coords.txt ]; then
        blastdbcmd -db ${db_path} -entry_batch extract_coords.txt -outfmt %f -out sampled_sequences.fasta
    else
        touch sampled_sequences.fasta
    fi

    # TODO encode IDs with seqkit https://bioinf.shenwei.me/seqkit/usage/#replace (Rename with number of record)

    """
}

workflow {
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
}