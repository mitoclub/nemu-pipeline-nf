# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2025-12-26

### Added
- Simultaneous processing pipeline for handling multiple protein sequences in parallel
- MSA (Multiple Sequence Alignment) pipeline with configurable modes: `auto_cdn`, `accurate_cdn`, `fast_cdn`, or `pure_mafft`
- Query sequence splitting process for batch processing
- Enhanced species filtering and BLAST result processing
- Configurable minimum sequence threshold (`min_seqs` parameter) before proceeding to MSA

### Improved
- BLAST step optimization with improved target sequence selection
- Species name parsing and species/relatives taxonomy resolution
- Sequence validation and unique sequence identification
- Filtering logic for sequence quality and diversity
- Overall pipeline beautification and code organization

### Fixed
- Species taxid parsing issues
- Outgroup-related bugs
- BLAST species and outgroup sequence encoding
- Multiple alignment pipeline flow

### Changed
- Refactored main pipeline to support simultaneous processing of multiple sequences
- Updated pipeline architecture for better maintainability

### Requirements
- Nextflow
- seqkit v2.9.0
- taxonkit v0.20.0
- BLAST+ 2.17.0
- Python 3
- mafft
- trimAl
- iqtree 2.2.0
- MACSE v2.07

## [1.0.0] - DATE TODO

### Initial Release
- Basic protein sequence analysis pipeline (NeMu)
- Single sequence file processing
- BLAST homology search
- Taxonomic data retrieval
- Phylogenetic analysis support 