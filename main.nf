#!/usr/bin/env nextflow

/*
 * Pipeline: Protein Sequence Analysis
 * Steps:
 * 1. Split Multi-FASTA & Assign Unique IDs
 * 2. Validate Protein sequences
 * 3. Parse Species Name
 * 4. Prepare Taxonomy (TaxonKit: Species & Relatives)
 * 5. TBLASTN & Filter (Single optimized step)
 */

// Global parameters
params.input        = "data/proteins.fa"
params.outdir       = "results"
params.db           = "/home/kpotoh/.nuc_db/dolphin"
params.taxdump      = "/home/kpotoh/.taxonkit" 
params.species_name = false  
params.gencode      = 2      
params.max_target_seqs = 2000 // Increased to ensure we capture enough hits for both species and outgroup

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
            rm family_all.taxid
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
 * 1. Merges TaxIDs (Species + Relatives)
 * 2. Runs Single TBLASTN
 * 3. Filters hits (Strict for Species, loose for Outgroup)
 * 4. Extracts sequences
 */
process TBLASTN_AND_FILTER {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'
    cpus 4 

    input:
    tuple val(id), path(query), path(species_txt), path(species_taxids), path(relatives_taxids)
    val db_path

    output:
    path "blast_result.tsv", optional: true
    path "sampled_sequences.fasta"
    path "filtering_log.txt"

    script:
    """
    # 1. Prepare Combined TaxID List
    cat ${species_taxids} ${relatives_taxids} > search.taxids
    
    # If list is empty (taxonkit failed), warn and exit
    if [ ! -s search.taxids ]; then
        echo "No TaxIDs found. Skipping BLAST." > filtering_log.txt
        touch sampled_sequences.fasta
        exit 0
    fi

    # 2. Run TBLASTN
    # We ask for 'staxids' to differentiate species from relatives later
    outfmt="6 saccver pident length qlen gapopen sstart send evalue bitscore sframe staxids"
    
    tblastn -query ${query} -db ${db_path} -db_gencode ${params.gencode} \
        -max_target_seqs ${params.max_target_seqs} \
        -evalue 0.0001 -num_threads ${task.cpus} \
        -taxidlist search.taxids \
        -outfmt "\$outfmt" \
        -out blast_result.tsv

    # 3. Python Script: Filter and Select Outgroup
    python3 -c "
import sys

def parse_taxids(filename):
    ids = set()
    try:
        with open(filename) as f:
            for line in f:
                if line.strip(): ids.add(line.strip())
    except FileNotFoundError:
        pass
    return ids

species_ids = parse_taxids('${species_taxids}')

species_hits = []
outgroup_hits = []

# Thresholds
MIN_IDENT_SPECIES = 80.0
MIN_COV_SPECIES = 0.5

MIN_IDENT_OUTGROUP = 60.0 # Looser for outgroup
MIN_COV_OUTGROUP = 0.5

try:
    with open('blast_result.tsv') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 11: continue
            
            # saccver pident length qlen gapopen sstart send evalue bitscore sframe staxids
            sacc = parts[0]
            pident = float(parts[1])
            length = float(parts[2])
            qlen = float(parts[3])
            bitscore = float(parts[8])
            staxids = parts[10].split(';') # Handle multiple taxids
            
            coverage = length / qlen
            
            # Check if ANY of the hit's taxids belong to our target species
            is_species = any(tid in species_ids for tid in staxids)
            
            hit_data = {
                'sacc': sacc,
                'sstart': parts[5],
                'send': parts[6],
                'sframe': parts[9],
                'bitscore': bitscore,
                'pident': pident,
                'line': line.strip()
            }

            if is_species:
                if pident >= MIN_IDENT_SPECIES and coverage >= MIN_COV_SPECIES:
                    species_hits.append(hit_data)
            else:
                # It is a relative (since we restricted blast to species+relatives)
                if pident >= MIN_IDENT_OUTGROUP and coverage >= MIN_COV_OUTGROUP:
                    outgroup_hits.append(hit_data)

except FileNotFoundError:
    print('No blast results found.')

# Sort Outgroup by bitscore (descending) to find the Closest Relative
outgroup_hits.sort(key=lambda x: x['bitscore'], reverse=True)

# Select Hits
final_coords = []
log_lines = []

log_lines.append(f'Found {len(species_hits)} valid species hits.')
for h in species_hits:
    # Format for blastdbcmd: ID start-end strand
    strand = 'minus' if int(h['sframe']) < 0 else 'plus'
    # Ensure start < end for blastdbcmd range extraction logic if needed, 
    # but blastdbcmd handles start>end as minus strand automatically if we format carefully.
    # However, legacy code logic: entry range strand
    # Let's use simple range format: start-end
    
    start, end = h['sstart'], h['send']
    final_coords.append(f\"{h['sacc']} {start}-{end} {strand}\")

log_lines.append(f'Found {len(outgroup_hits)} potential outgroups.')
if outgroup_hits:
    best_out = outgroup_hits[0]
    log_lines.append(f\"Selected Outgroup: {best_out['sacc']} (Ident: {best_out['pident']}%, Score: {best_out['bitscore']})\")
    
    strand = 'minus' if int(best_out['sframe']) < 0 else 'plus'
    start, end = best_out['sstart'], best_out['send']
    final_coords.append(f\"{best_out['sacc']} {start}-{end} {strand}\")
else:
    log_lines.append('WARNING: No valid outgroup found.')

with open('extract_coords.txt', 'w') as f:
    for line in final_coords:
        f.write(line + '\\n')

with open('filtering_log.txt', 'w') as f:
    f.write('\\n'.join(log_lines) + '\\n')
    "

    # 4. Extract Sequences using blastdbcmd
    if [ -s extract_coords.txt ]; then
        # blastdbcmd -entry_batch expects: ID range strand? 
        # Actually standard input for -entry_batch is just ID, or ID range. 
        # To handle strand, it's safer to use the 'awk' trick from legacy or handle extraction carefully.
        # But blastdbcmd -entry_batch supports "seq_id range strand" in newer versions? 
        # Standard: id start-end strand
        
        blastdbcmd -db ${db_path} -entry_batch extract_coords.txt -outfmt %f -out sampled_sequences.fasta
    else
        touch sampled_sequences.fasta
    fi
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
            if (seq.length() > 0 && (dna_count / seq.length()) > 0.95) return false
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
    SAVE_QUERY(PREPARE_TAXONOMY.out)
    
    TBLASTN_AND_FILTER(SAVE_QUERY.out, params.db)
}