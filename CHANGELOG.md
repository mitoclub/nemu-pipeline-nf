# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.3] - 2026-08-14

### Fixed

- Outgroup remapping in `ENCODE_AND_RMDUP` now matches the original FASTA header/ID instead of always relabeling the last sequence (broke nucleotide input when OUTGRP was not last, and could rename `seq_1` plus `seq_10+` via unanchored regex)
- Nucleotide outgroup detection now inspects FASTA headers rather than a raw substring of the whole file
- `seqkit rmdup` no longer drops OUTGRP when its sequence is identical to an ingroup record
- `FILTER_AND_EXPORT` always sets `OUTGRP_ID` so missing outgroups continue instead of failing the process
- Mutation-count filter after `MUT_EXTRACTION` uses the number of table rows, not file size in bytes (`--min-muts`)
- Nucleotide `--input` globs that match nothing now fail (`checkIfExists`)
- Coding MSA gap cleaning now drops codon columns (multiples of 3) so the reading frame is preserved
- N-content warning uses `seqkit fx2tab --base-content N` instead of GC/quality columns
- Spectra calculation passes `--proba` when mutation extraction used probabilities
- Sequences are uppercased before encoding; BLAST DB existence accepts `.ndb`, `.nin`, or `.nal`
- Example E. coli FASTA outgroup headers renamed from `>outgroup` to `>OUTGRP` to match the documented default

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

## [1.0.0] - 2026-01-28

### Initial Release

- Single input file processing
- Singularity container is the main environment
- Core functionality described in the NeMu paper
