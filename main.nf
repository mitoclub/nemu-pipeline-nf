#!/usr/bin/env nextflow

/*
 * Pipeline: Protein Sequence Analysis
 * Steps:
 * 1. Split Multi-FASTA & Assign Unique IDs
 * 2. Validate Protein sequences
 * 3. Parse Species Name
 * 4. Prepare Taxonomy (TaxonKit: Species & Relatives)
 * 5. Save Query & Metadata
 * 6. Run TBLASTN (Restricted by TaxID)
 */

// Global parameters
params.input        = "data/proteins.fa"
params.outdir       = "results"
params.db           = "/home/kpotoh/.nuc_db/dolphin_db"  // Path to BLAST database
params.taxdump      = "/home/kpotoh/.taxonkit" // Directory containing nodes.dmp and names.dmp
params.species_name = false  // Optional: Provide species name directly
params.gencode      = 2      // Genetic code for translation
params.max_target_seqs = 1000  // Max target sequences for BLAST

log.info """\
    P R O T E I N   P I P E L I N E
    ===================================
    input file   : ${params.input}
    outdir       : ${params.outdir}
    blast db     : ${params.db}
    taxdump dir  : ${params.taxdump}
    species_name : ${params.species_name ?: 'Auto-detect from headers'}
    """
    .stripIndent()

/*
 * PROCESS: Parse Species Name
 * Uses input header to determine species.
 */
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

/*
 * PROCESS: Prepare Taxonomy
 * Uses taxonkit to derive TaxIDs for the species and its relatives.
 */
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

    # Clean species name (ensure spaces for taxonkit)
    CLEAN_NAME=\$(echo "${species_name}" | tr '_' ' ')

    echo "Deriving TaxIDs for: \$CLEAN_NAME"

    # 1. Get TaxID
    SPEC_ID=\$(echo "\$CLEAN_NAME" | taxonkit name2taxid | cut -f2)

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
            # 4. Get Family members and exclude self to find relatives
            taxonkit list --ids \$FAMILY_ID --indent "" | head -n -1 > family_all.taxid
            grep -vFf species.taxid family_all.taxid > relatives.taxid
            # rm family_all.taxid
        fi
    fi
    """
}

/*
 * PROCESS: Save Query & Metadata
 */
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
 * PROCESS: TBLASTN
 * Uses available CPUs to speed up search.
 */
process TBLASTN {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'
    
    // Default to 4 CPUs per BLAST task. 
    // Nextflow will run as many tasks in parallel as your machine permits (Total CPUs / 4).
    cpus 4 

    input:
    tuple val(id), path(query), path(species_txt), path(species_taxid_list), path(relatives_taxid_list)
    val db_path

    output:
    path "blast_output_species.tsv"

    script:
    """
    outfmt="6 saccver pident length qlen gapopen sstart send evalue bitscore sframe staxids"

    ARGS=""
    if [ -s "${species_taxid_list}" ]; then
        ARGS="-taxidlist ${species_taxid_list}"
        echo "Running with taxid restriction..."
    fi

	tblastn -query ${query} -db ${db_path} -db_gencode ${params.gencode} \
        -max_target_seqs ${params.max_target_seqs} \
		-evalue 0.00001 -num_threads ${task.cpus} \$ARGS \
		-outfmt "\$outfmt" -out blast_output_species.tsv
    """
}

workflow {
    // Initialize a counter for unique IDs
    def seq_counter = 0

    // 1. Prepare Channels
    raw_sequences = Channel.fromPath(params.input)
        .splitFasta(record: [id: true, header: true, seqString: true])
        .filter { record ->
            // Validation Logic
            def seq = record.seqString.toUpperCase()
            if (seq =~ /[EFILPQZ]/) return true
            def dna_count = seq.count('A') + seq.count('C') + seq.count('G') + seq.count('T') + seq.count('N')
            if (seq.length() > 0 && (dna_count / seq.length()) > 0.95) return false
            return true
        }
        .map { record ->
            // Increment counter for every sequence passed
            def count = ++seq_counter
            def clean_original_id = record.id.split()[0].replaceAll(/[^a-zA-Z0-9\.]/, '_')            
            def unique_id = "${count}__${clean_original_id}"
            [unique_id, record.header, record.seqString]
        }

    // 2. Parse Species Name
    PARSE_SPECIES_NAME(raw_sequences, params.species_name)

    // 3. Prepare Taxonomy
    PREPARE_TAXONOMY(PARSE_SPECIES_NAME.out, params.taxdump)

    // 4. Save Query Files
    SAVE_QUERY(PREPARE_TAXONOMY.out)

    // 5. Run TBLASTN
    TBLASTN(SAVE_QUERY.out, params.db)
}