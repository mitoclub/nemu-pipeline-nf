# NeMu-pipeline

A Nextflow-based bioinformatics pipeline for sequence analysis and phylogenetic inference.

<!-- example repo - https://github.com/cbcrg/unistrap/tree/master -->

## Workflow schematic representation

Image TODO


## Dependencies

All dependencies are specified in [environment.yml](./environment.yml) and include:

- **Nextflow** (with **Java**)
- **Python**
- **IQtree**
- **SeqKit**
- **TaxonKit**
- **BLAST+**
- **MAFFT** and **MACSE**

## Installation

### Conda

Install all dependencies using conda or mamba. Conda and mamba can be installed following [these instructions](https://github.com/conda-forge/miniforge)

```bash
conda create -f environment.yml -n nemu-pipeline --yes
conda activate nemu-pipeline
```

### Docker

Build the container to run the pipeline inside it. You can either run the entire pipeline with Nextflow inside the container or use the container as an isolated environment for the Nextflow pipeline (requires manually installing Nextflow 25.10 with Java OpenJDK 17).

```bash
docker build -t nemu-pipeline:latest .
#docker run -it nemu-pipeline:latest nextflow run /app/main.nf ...
```

## Usage

Clone the repository and run the pipeline with your input sequences:

```bash
git clone <repository-url>
cd nemu-pipeline

nextflow run main.nf \
  -process.cpus=20 \
  -resume \
  -with-trace \
  --input "sample_input_ecoli_head/*.fasta" \
  --input-type nucleotide_coding \
  -o results_ecoli \
  --gencode 1 \
  --outgroup-id outgroup
```

**Parameters:**
- `--input`: Path pattern to input FASTA files
- `--input-type`: Type of input sequence (e.g., nucleotide_coding)
- `-o`: Output directory for results
- `--gencode`: Genetic code table to use
- `--outgroup-id`: Identifier for the outgroup sequence


## Command line options

TODO


## TODO

- [ ] Add parsing of taxon IDs from FASTA headers (configurable option)
- [ ] Add comprehensive tests for pipeline (Nextflow-based test files)
- [ ] Run and validate tests
- [ ] Finalize pipeline configuration
- [ ] Prepare and optimize container image


## Testing

Run the test script from the repository directory to validate the pipeline:

```bash
bash test.sh
```

