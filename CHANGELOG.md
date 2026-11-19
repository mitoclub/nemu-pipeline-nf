# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2026-11-19

### Added

BLAST identity and coverage, alignment gap fractions, the ancestral-state probability cutoff, the TreeShrink quantile, the 192-class observation minimum, and the outgroup branch rank are parameters. The defaults are the previous fixed values: `--minIdentSpecies 80`, `--minCovSpecies 0.5`, `--minIdentOutgroup 70`, `--minCovOutgroup 0.5`, `--maxGapSeq 0.5`, `--maxGapSite 0.5`, `--probaCutoff 0.3`, `--treeShrinkQuantile 0.1`, `--mnum192 16`, and `--outgroupTopN 5`. Identity values are percents. Coverage, gap fractions, the probability cutoff, and the TreeShrink quantile are fractions from 0 to 1.

### Changed

`--inputType` accepts only `protein`, `cds`, and `noncds`. Older names (`prot`, `CDS`, `nucleotide_coding`, `NONCDS`, `nucleotide_noncoding`) stop the run and name the value to use instead. `noncds` still forces MAFFT.

### Fixed

A nucleotide run no longer stops at the start because the taxonomy output is empty. `parsed_taxonomy.txt` is published only for protein input.

## [1.2.0] - 2026-09-07

### Added

`--minMuts` (default 5) sets how many reconstructed mutations a gene needs before a spectrum is derived. The old check compared the observation file’s size in bytes, so real spectra with few mutations were dropped and a large file could pass with almost no mutations.

### Fixed

Outgroup, alignment, and spectrum steps now follow the sequence the run actually claims to use. Nucleotide outgroups are taken from FASTA headers and renamed from that record, including when OUTGRP is not the last sequence and when encoded ids such as `seq_1` and `seq_10` would previously be rewritten together. A missing outgroup is reported and the protein path continues. Coding alignments drop gappy columns in whole codons so the reading frame stays intact, a short outgroup is kept through the length filter, and sequences are uppercased before encoding. A nucleotide `--input` that matches no files fails immediately. The high-N warning reads the N percentage, spectrum calculation receives `--proba` when mutation extraction used probabilities, and a BLAST database is accepted when `.ndb`, `.nin`, or `.nal` is present.

## [1.1.0] - 2026-03-11

### Added

- Query sequence splitting for batch processing
- Simultaneous processing pipeline for handling multiple inputs in parallel
- MSA (Multiple Sequence Alignment) pipeline with configurable modes: `auto`, `macse`, `mafft_macse`, or `mafft`
- Configurable minimum sequence threshold (`minSeqs` parameter) before proceeding to MSA
- Tree visualization process using nw_display
- Possibility to run without outgroup

### Improved

- Taxonomy ID retrieval with enhanced error handling
- BLAST step optimization with improved target sequence selection
- Filtering logic for sequence quality and relevance with nextflow functions
- Enhanced species filtering and BLAST result processing
- Separated workflows: nemuCore and blastHead for better modularity
- Overall pipeline beautification and code organization

### Fixed

- Species taxid parsing issues (taxonkit integration)
- Outgroup-related bugs
- BLAST species and outgroup sequence encoding
- Channel output wiring issue in nucleotide input path (`CHECK_INPUT_TYPE`) that caused runtime aborts
- Outgroup remapping logic in `ENCODE_AND_RMDUP` to avoid incorrect fallback to the last sequence
- Sequence-threshold checks to use `>= minSeqs` consistently

## [1.0.1] - 2024-05-08

### Initial Release

- Single input file processing
- Singularity container is the main environment
- Core functionality described in the NeMu paper
