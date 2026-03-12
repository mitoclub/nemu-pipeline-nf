# NeMu-pipeline

A Nextflow-based bioinformatics pipeline for mutational spectra reconstruction based on sequences using phylogenetic methods.

<!-- example repo - https://github.com/cbcrg/unistrap/tree/master -->

## Workflow schematic representation

![scheme](./docs/scheme.jpg)

## Features

- Authomatic sequences retrieval using tblastn and nucleotide database (when the input is protein sequence)
- Outlier sequences removal 
- Mutation polarization using outrgoup based phylogenetic tree rooting
- Accounting of mutations probabilities according to ancestral states probabilities
- Sampling of mutations along the tree for variance estimation
- Spectra derivation for all synonymous and synonymous fourfold (syn4f) sites separately
- Spectra derivation according to site rate category (optional)

## Dependencies

All dependencies are specified in [environment.yml](./environment.yml) and include:

- **Nextflow** (with **Java**)
- **BLAST+**
- **MAFFT** and **MACSE**
- **SeqKit** and **TaxonKit**
- **IQtree**
- **TreeShrink**
- **Python** (**PyMutSpec** package)

## Installation

### Conda

Clone the repository. Install the dependencies using conda or mamba. 

Conda and mamba can be installed following [these instructions](https://github.com/conda-forge/miniforge).

```bash
git clone https://github.com/mitoclub/nemu-pipeline-nf.git
cd nemu-pipeline

conda env create -f environment.yml -n nemu --yes
# OR
# mamba create -f environment.yml -n nemu --yes

conda activate nemu-pipeline
```

### Docker

Build the container to run the pipeline inside it. You can either run the entire pipeline with Nextflow inside the container or use the container as an isolated environment for the Nextflow pipeline (requires manually installing Nextflow 25.10 with Java OpenJDK 17).

```bash
docker build -t nemu-pipeline:latest .
#docker run -it nemu-pipeline:latest nextflow run /app/main.nf ...
```

## Usage

### When input is nucleotide multi-fasta file

Sequences must be orthologous and may be aligned (--aligned parameter).

It's possible to pass several input fasta files using "*". Note that it's required to use caveats for the filename.

```bash
nextflow run main.nf \
  -o results_ecoli \
  -process.cpus=20 \
  --input-type nucleotide_coding \
  --input "sample_input_ecoli/*.fasta" \
  --gencode 1 \
  --outgroup-id outgroup
```

### When input is protein sequences

```bash
nextflow run main.nf \
  -o results_test \
  -process.cpus=20 \
  -resume \
  -with-trace \
  --input-type protein \
  --input "test_data/test_proteins.fasta" \
  --gencode 2 \
  --db path_to_nuc_database \
  --taxdump "$HOME/.taxonkit"
```

Recomendations for comparative-species analysis: --branch-spectra, --model, OUTGRP, consCatCutoff etc. TODO

## Command line options

It's possible to use command line parameters or change nextflow.config file.

```txt
Usage: nextflow run main.nf --input <input_fasta> [options]

TODO update to latest

Main options:
    --input FILE            Input FASTA file (required)
                            If input type is protein, a multi-FASTA with one or several sequences 
                            required (header format: ">ID [Species name]"). If input type 
                            is nucleotide, one or several fasta files with orthologous 
                            sequences (including outgroup) required
    --input_type STRING     Type of input sequences: 
                            protein, nucleotide_coding, nucleotide_noncoding (default: nucleotide_coding)
    --gencode NUM           Genetic code table (default: 1)
                            Used for codon-aware alignment and annotation of mutations

Required options for nucleotide input:
    --outgroup-id STRING    Outgroup sequence ID (default: OUTGRP)
                            Specify outgroup sequence ID in the alignment for rooting the tree
    --aligned BOOL          Input sequences are pre-aligned (default: false)

Required options for protein input:
    --db PATH               BLAST database path (required for protein input)
    --taxdump DIR           Taxdump directory path (required for protein input)
                            TODO integrate to the container
    --species-name STRING   Override species name. Useful when you work with proteins 
                            from single species
    --max-target-seqs NUM   Max target sequences for BLAST (default: 2000)
                            tblastn will collect no more than this number of sequences

Nextflow options:
    -with-report FILE       Generate execution report
    -with-trace FILE        Generate execution trace
    -with-timeline FILE     Generate execution timeline
    -output-dir DIR         Specify output directory (default: ./results)

Common options:
    --threads NUM           Number of threads to use (default: 1) TODO delete if not needed
    --save-intermeds BOOL   Save intermediate files TODO implement
    --help                  Show this help message and exit

Options for MSA & Phylogeny:
    --msa-mode STRING       MSA mode: auto, macse, mafft_macse, mafft (default: auto)
    --min-seqs NUM          Minimum number of sequences to proceed phylogenetic inference (default: 4)
    --treefile FILE         Input tree file (optional; default: build tree de novo)
    --model STRING          IQ-TREE substitution model (default: GTR+FO+G4+I)
    --model-asr STRING      ASR substitution model (default: GTR+FO+G4+I)
    --run-treeshrink BOOL   Run TreeShrink to prune long branches (default: true)

Options for Mutation Extraction:
    --cons-cat-cutoff NUM   Conservation category cutoff for mutation extraction (default: 0)
                            0 = no cutoff; only mutations in sites with rate category 
                            less than or equal to this value will be used
    --proba-arg BOOL        Use probabilities in mutation extraction (default: true)
    --uncertainty-coef BOOL Use phylogeny uncertainty coefficient in mutation extraction 
                            (default: true)

Options for Mutation Spectra Derivation:
    --plot BOOL             Generate barcharts of mutation spectra (default: true)
    --internal BOOL         Derive spectra for internal branches (default: false)
    --terminal BOOL         Derive spectra for terminal branches (default: false)
    --branch-spectra BOOL   Derive spectra for individual branches (default: false)
```

## TODO

- [ ] Add parsing of taxon IDs from FASTA headers (configurable option)
- [ ] Add comprehensive tests for pipeline (Nextflow-based test files)
- [ ] Run and validate tests
- [ ] Finalize pipeline configuration
- [ ] Prepare and optimize container image

<!-- 
## Testing

Run the test script from the repository directory to validate the pipeline:

```bash
bash test.sh
```
 -->
