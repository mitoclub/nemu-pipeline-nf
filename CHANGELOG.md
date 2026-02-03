# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2026-02-03

### Added

- Query sequence splitting for batch processing
- Simultaneous processing pipeline for handling multiple inputs in parallel
- MSA (Multiple Sequence Alignment) pipeline with configurable modes: `auto_cdn`, `accurate_cdn`, `fast_cdn`, or `pure_mafft`
- Configurable minimum sequence threshold (`min_seqs` parameter) before proceeding to MSA
- Profiles in the config file for different environments (e.g. conda, singularity, docker)
- Tree visualization process using nw_display

### Improved

- Taxonomy ID retrieval with enhanced error handling
- BLAST step optimization with improved target sequence selection
- Filtering logic for sequence quality and relevance with nextflow functions
- Enhanced species filtering and BLAST result processing
- Separated workflows: nemuCore and blastHead for better modularity
- Overall pipeline beautification and code organization

### Fixed

- Species taxid parsing issues (taxonkit integration)
- Outgroup-related bugs and possibility to run without outgroup
- BLAST species and outgroup sequence encoding

## [1.0.0] - DATE TODO

### Initial Release

- Single input file processing
- Singularity container is the main environment
- Functionality described in the paper TODO
