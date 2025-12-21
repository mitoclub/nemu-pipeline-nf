#!/usr/bin/env nextflow

/*
 * Pipeline: Protein Sequence Analysis
 * Steps:
 * 1. Split Multi-FASTA
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

    m_os = re.search(r'OS=([a-zA-Z0-9_ ]+)', header)
    if m_os:
        species = m_os.group(1).strip()
    else:
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
    # Configure taxonkit to use the provided dump files
    export TAXONKIT_DB=${taxdump_dir}

    # Clean species name (replace underscores with spaces for taxonkit lookup if needed)
    # Most taxonkit lookups work better with spaces: "Homo_sapiens" -> "Homo sapiens"
    CLEAN_NAME=\$(echo "${species_name}" | tr '_' ' ')

    echo "Deriving TaxIDs for: \$CLEAN_NAME"

    # 1. Get TaxID for the Species
    # name2taxid returns: "Name <tab> TaxID"
    SPEC_ID=\$(echo "\$CLEAN_NAME" | taxonkit name2taxid | cut -f2)

    if [ -z "\$SPEC_ID" ]; then
        echo "WARNING: TaxID not found for \$CLEAN_NAME"
        touch species.taxid relatives.taxid
    else
        # 2. Get all downstream TaxIDs for this Species (Subspecies, strains)
        # This file is used to restrict TBLASTN search
        taxonkit list --ids \$SPEC_ID --indent "" > species.taxid

        # 3. Find the Family ID to determine relatives
        # reformat -F -f "{k}\t{t}" outputs rank and taxid for the lineage
        FAMILY_ID=\$(echo \$SPEC_ID | taxonkit lineage | taxonkit reformat -F -f "{k}\t{t}" | grep -P "^family\t" | cut -f2)

        if [ -z "\$FAMILY_ID" ]; then
            echo "WARNING: Family rank not found for \$CLEAN_NAME (ID: \$SPEC_ID)"
            touch relatives.taxid
        else
            # 4. Get all members of the Family
            taxonkit list --ids \$FAMILY_ID --indent "" > family_all.taxid
            
            # 5. Exclude the Species IDs from the Family IDs to get strictly relatives (Outgroup candidates)
            # grep -vFf uses species.taxid as a pattern file to exclude lines from family_all.taxid
            grep -vFf species.taxid family_all.taxid > relatives.taxid
            
            # Cleanup temp file
            rm family_all.taxid
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
 * Uses -taxidlist to speed up search if the list is available.
 */
process TBLASTN {
    tag "$id"
    publishDir "${params.outdir}/${id}", mode: 'copy'

    input:
    tuple val(id), path(query), path(species_txt), path(species_taxid_list), path(relatives_taxid_list)
    val db_path

    output:
    path "blast_result.txt"

    script:
    """
    # Check if DB exists
    if ! ls ${db_path}* 1> /dev/null 2>&1; then
        echo "WARNING: BLAST DB not found" > blast_result.txt
        exit 0
    fi

    ARGS="-query ${query} -db ${db_path} -outfmt 6"

    # If species.taxid exists and is not empty, use it to restrict search
    if [ -s "${species_taxid_list}" ]; then
        ARGS="\$ARGS -taxidlist ${species_taxid_list}"
        echo "Running TBLASTN with taxid restriction..."
    else
        echo "Running TBLASTN without taxid restriction (list empty or missing)..."
    fi

    tblastn \$ARGS -out blast_result.txt
    """
}

workflow {
    def seq_counter = 0

    // 1. Prepare Channels
    raw_sequences = Channel.fromPath(params.input)
        .splitFasta(record: [id: true, header: true, seqString: true])
        .filter { record ->
            // Simple validation
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

    // 2. Parse Species Name
    PARSE_SPECIES_NAME(raw_sequences, params.species_name)

    // 3. Prepare Taxonomy (Generates taxid lists)
    PREPARE_TAXONOMY(PARSE_SPECIES_NAME.out, params.taxdump)

    // 4. Save Query Files
    SAVE_QUERY(PREPARE_TAXONOMY.out)

    // 5. Run TBLASTN
    TBLASTN(SAVE_QUERY.out, params.db)
}