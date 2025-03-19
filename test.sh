#!/bin/bash

# Custom sequences
seq_nucs=../example_data/CO1__Stenella_longirostris.fna
if [ -f nemu-core.nf ] && [ -f test-core.config ] && [ -f $seq_nucs ]; then
    nextflow -c test-core.config run -resume nemu-core.nf --sequence $seq_nucs --outdir test_output1
else
    echo you must be in the pipeline directory to access nemu-core.nf and test-core.config files
    exit 1
fi

# # MIDORI2 database
# seq_protein=../example_data/CYTB__Stenella_longirostris.faa
# if [ -f nemu.nf ] && [ -f test.config ] && [ -f $seq_protein ]; then
#     nextflow -c test.config run -resume nemu.nf --sequence $seq_protein --outdir test_output2
# else
#     echo you must be in the pipeline directory to access nemu.nf and test.config files
#     exit 1
# fi

# # NT database
# seq_protein=../example_data/CYTB__Stenella_longirostris.faa
# db=/db/nt
# genus_taxid=9734
# if [ -f nemu.nf ] && [ -f test.config ] && [ -f $seq_protein ]; then
#     nextflow -c test.config run -resume nemu.nf --sequence $seq_protein --outdir test_output3 \
#     --DB $db --genus_taxid $genus_taxid
# else
#     echo you must be in the pipeline directory to access nemu.nf and test.config files
#     exit 1
# fi

# # NT database error
# seq_protein=../example_data/CYTB__Stenella_longirostris.faa
# db=/db/nt
# if [ -f nemu.nf ] && [ -f test.config ] && [ -f $seq_protein ]; then
#     nextflow -c test.config run -resume nemu.nf --sequence $seq_protein --outdir test_output3 \
#     --DB $db
# else
#     echo you must be in the pipeline directory to access nemu.nf and test.config files
#     exit 1
# fi

echo DONE
