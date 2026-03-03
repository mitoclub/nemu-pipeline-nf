# The pipeline for neutral mutation spectra evaluation based on evolutionary data

# USe this https://github.com/cbcrg/unistrap/tree/master to 

## Quick start

```bash
git clone TODO
cd nemu-pipeline

./main.nf --input "sample_input_ecoli_head/*.fasta" --input-type nucleotide_coding -o results_ecoli --gencode 1 --outgroup-id outgroup

```


## TODO

- add parsing of taxids from fasta headers (using additional option or not) (NEED TO TEST NOW)
- add tests for nemu-core pipeline (without bash script, just nextflow test files)
- run tests
- fix config
- prepare container

## 2 pipeline versions

1. [NeMu pipeline including tblastn-head](./nemu.nf) - input is single protein sequence that will be used by tblastn to search homologous nucleotide sequences in selected database. After this step NeMu-core executes with phylogeny and spectra inference
2. [NeMu-core pipeline](./nemu-core.nf) - input is multifasta of nucleotide sequences, that used for phylogeny and spectra inference

## Config examples

1. [Config for comparative species analysis](./comp_sp.config) - on many species
2. [Config for intraspecies analysis](./single_sp.config) - on single species

    - Don't forget to change process.container and singularity.runOptions parameters according to execution environment

## Test

TODO

Run the script (from this directory) to test the pipeline on your computer. Don't forget to change path to singularity container and runOptions in the config file.

**Important about runOptions:** if you run the pipeline from disk that don't contain your $HOME directory, you must write in the runOptions `--bind $MOUNT_PATH`, where $MOUNT_PATH is the mount point of the disk from that you want to execute the pipeline (`--bind /scratch` in my case). If you run the pipeline from any subdirectory of your $HOME, delete this runOptions from config file.

```bash
bash test.sh
```


## Build and activate environment

### Mamba

```bash
mamba create -f environment.yml -n nemu-pipeline --yes
mamba activate nemu-pipeline
```

### Apptainer/Singularity
```bash
apptainer build nemu-pipeline.sif nemu-pipeline.def
apptainer exec nemu-pipeline.sif nextflow run /app/main.nf ...
```

### Docker

```bash
docker build -t nemu-pipeline:latest .
docker run -it nemu-pipeline:latest nextflow run /app/main.nf ...
```