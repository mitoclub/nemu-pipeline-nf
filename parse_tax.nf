#!/usr/bin/env nextflow

params.input            = ""
params.taxdump          = ""
params.speciesName      = ""                            // Override species name

process PREPARE_TAXONOMY {
    tag "$id"

    input:
    tuple val(id), val(species_name), val(sequence)
    path taxdump_dir

    output:
    tuple val(id), path("query.fa"), path("species.taxid"), path("relatives.taxid"), val(species_name)

    script:
    """
    if [ "${species_name}" = "unknown_species" ]; then
        touch species.taxid relatives.taxid query.fa
        exit 0
    fi

    export TAXONKIT_DB=${taxdump_dir}
    CLEAN_NAME=\$(echo "${species_name}" | tr '_' ' ')

    if [[ \$CLEAN_NAME =~ ^[0-9]+\$ ]]; then
        echo "Assuming provided species name is TaxID: \$CLEAN_NAME"
        SPEC_ID=\$CLEAN_NAME
        echo "Check existance of TaxID \$SPEC_ID"
        echo \$SPEC_ID | taxonkit lineage -c > given_taxid_lineage.txt
        status_code=\$(cut -f2 given_taxid_lineage.txt)
        if [ \$status_code != \$SPEC_ID ]; then
            echo "WARNING: Provided TaxID \$SPEC_ID not found in taxonomy database"
            touch species.taxid relatives.taxid query.fa
            exit 0
        fi
    else
        echo "Deriving TaxIDs using taxonkit for: \$CLEAN_NAME"
        
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
    fi
    
    # 2. Get downstream TaxIDs
    taxonkit list --ids \$SPEC_ID --indent "" | head -n -1 > species.taxid

    # 3. Find Family ID  TODO do this single time like in L117 (collect species,genus,family)
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
    echo ">${id}" > query.fa
    echo "${sequence}" >> query.fa
    """
}

process PREPARE_TAXONOMY_NEW {
    input:
    path species_list
    path taxdump_dir

    output:
    path "taxonomy_final.txt"

    script:
    """
    export TAXONKIT_DB=${taxdump_dir}

    # Split numeric taxids and species names in one pass.
    awk '
        /^[0-9]+\$/ { print > "taxids.txt"; next }
        NF { gsub(/_/, " "); print > "species_names.txt" }
    ' ${species_list}
    [ -f taxids.txt ] || : > taxids.txt
    [ -f species_names.txt ] || : > species_names.txt

    # Derive taxids only when species names are present.
    if [ -s species_names.txt ]; then
        taxonkit name2taxid --show-rank species_names.txt > species_info.txt || true
    else
        : > species_info.txt
    fi

    # Merge taxids from direct numeric inputs and resolved species names.
    cut -f2 species_info.txt > taxids_from_names.txt
    cat taxids.txt taxids_from_names.txt | awk '/^[0-9]+\$/' | sort -n -u > all_taxids.txt

    if [ ! -s all_taxids.txt ]; then
        : > taxonomy_lineages.txt
        echo '{}' > taxonlist.json
        : > taxonomy_final.txt
        exit 0
    fi

    taxonkit lineage all_taxids.txt | taxonkit reformat -t -f "{f},{s}" > taxonomy_lineages.txt

    # Prepare a single id list for taxonkit list --ids.
    tx_list_to_parse=\$(cut -f4 taxonomy_lineages.txt | tr ',' '\n' | awk '/^[0-9]+\$/' | sort -nu | paste -sd ',' -)
    if [ -z "\$tx_list_to_parse" ]; then
        tx_list_to_parse=\$(paste -sd ',' all_taxids.txt)
    fi

    if [ -n "\$tx_list_to_parse" ]; then
        taxonkit list --ids "\$tx_list_to_parse" --json > taxonlist.json
    else
        echo '{}' > taxonlist.json
    fi

    # extract paths to leaf nodes (species) for lineage parsing in the next step
    cat taxonlist.json | jq -r 'paths(objects | select(length == 0)) | join(",")' > paths.csv

    cut -f4 taxonomy_lineages.txt > fam_sp_taxids.csv

    : > descendants.txt

    # iterate over family and species taxids and extract their lineages from taxonlist.json
    while IFS=, read -r fam_taxid sp_taxid; do
        echo "Processing \$fam_taxid and \$sp_taxid"
        echo \$sp_taxid > sp_lineage_\${sp_taxid}.txt
        cat taxonlist.json | jq --arg tax \$sp_taxid -r '.[\$tax]  | paths(objects | select(length == 0)) | join("\n")' | sort -nu >> sp_lineage_\${sp_taxid}.txt
        cat taxonlist.json | jq --arg tax \$fam_taxid -r '.[\$tax] | paths(objects | select(length == 0)) | join("\n")' | sort -nu > fam_lineage_\${fam_taxid}.txt
        
        grep -vf sp_lineage_\${sp_taxid}.txt fam_lineage_\${fam_taxid}.txt > fam_lineage_excl_sp_\${fam_taxid}.txt

        sp_lineage_lst=\$(paste -sd ',' sp_lineage_\${sp_taxid}.txt)
        fam_lineage_lst=\$(paste -sd ',' fam_lineage_excl_sp_\${fam_taxid}.txt)

        paste <(echo "\$sp_lineage_lst") <(echo "\$fam_lineage_lst") >> descendants.txt

    done < fam_sp_taxids.csv
    
    # TODO name in the input can differ from name in the taxonomy due to db versions
    # need to explicitly add name from the input to the firsh column
    paste <(cat taxonomy_lineages.txt) <(cat descendants.txt) >> taxonomy_final.txt
    """
}

process COMPARE_LINEAGES {
    tag "$id"
    errorStrategy 'ignore'

    input:
    tuple val(id), path(query), path(species_taxids), path(relatives_taxids), val(species_name)
    path taxonomy_final

    output:
    tuple val(id), path("lineage_comparison.txt")

    script:
    """
    # extract species and family taxids for the given species name from taxonomy_final
    awk -F'\t' -v name="$species_name" '
        \$1 == name || \$2 ~ name {print}
    ' $taxonomy_final > matched_species.txt
    
    fam_taxid=\$(cut -f4 matched_species.txt | cut -f1 -d',')

    cut -f5 matched_species.txt | tr "," "\n" | sort -n > new_sp_taxids.txt
    cut -f6 matched_species.txt | tr "," "\n" > new_fam_taxids.txt

    # prepare old taxids for comparison
    sort -nu $species_taxids > old_sp_taxids.txt
    sort -nu $relatives_taxids | grep -v "\$fam_taxid" > old_fam_taxids.txt
    
    touch lineage_comparison.txt
    diff -q old_sp_taxids.txt new_sp_taxids.txt > lineage_comparison.txt | true
    diff -q old_fam_taxids.txt new_fam_taxids.txt >> lineage_comparison.txt | true
    """
}

workflow {
    main:
    input_fasta = params.input
    taxdump = params.taxdump
    species_name = params.speciesName


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
        .unique().flatten().collectFile(name: 'sample.txt', newLine: true, sort: true)
    
    if (params.verbose) {
    species_lst.subscribe { file ->
            println "Parsed species names are saved to file: $file\n\n"
            // println "File content is:\n${file.text}"
        }
    }
            
    PREPARE_TAXONOMY_NEW(species_lst, taxdump)
    taxonomy_final = PREPARE_TAXONOMY_NEW.out.first()

    PREPARE_TAXONOMY(raw_sequences, taxdump)
        // subscribe { file ->
        //     println "Lineage comparison output saved to: $file"
        //     println "File content is:\n${file.text}"
        // }

    COMPARE_LINEAGES(PREPARE_TAXONOMY.out, taxonomy_final)

    // iterate over lineage comparison outputs and print warnings if there are differences
    COMPARE_LINEAGES.out.subscribe { tuple ->
        def (id, comparison_file) = tuple
        if (comparison_file.text.trim()) {
            log.warn "Lineage comparison for ${id} shows differences:\n${comparison_file.text}"
        }
    }

    // publish:
    // compares = COMPARE_LINEAGES.out
}

// output {
//     compares {
//         path { sample -> "${sample[0]}/" }
//     }
// }
