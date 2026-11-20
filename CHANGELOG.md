# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2026-11-19

### Added

- `--minIdentSpecies` (default `80`) and `--minCovSpecies` (default `0.5`) set the percent identity and query coverage required for a species BLAST hit. These values were previously fixed.
- `--minIdentOutgroup` (default `70`) and `--minCovOutgroup` (default `0.5`) set the same cutoffs for the outgroup BLAST hit.
- `--maxGapSeq` and `--maxGapSite` (default `0.5`) drop sequences or alignment columns whose gap fraction is above that value.
- `--minProba` (default `0.3`) sets the minimum ancestral-state probability. Use `--proba false` to ignore probabilities.
- `--treeShrinkQuantile` (default `0.1`) sets the TreeShrink branch-length quantile.
- `--minMuts192` (default `16`) sets the minimum observations required for each 192-class context.
- `--outgroupTopN` (default `5`) keeps the outgroup when its branch is among this many longest branches.
- Identity values are percents. Coverage, gap fractions, the probability cutoff, and the TreeShrink quantile are fractions from 0 to 1.

### Changed

- `--inputType` accepts only `protein`, `cds`, and `noncds`.
- Older names (`prot`, `CDS`, `nucleotide_coding`, `NONCDS`, `nucleotide_noncoding`) stop the run and name the value to use instead.
- `noncds` still forces MAFFT.
- Cryptic parameter names from 1.x are renamed: `--maxTargetSeqs` to `--maxHits`, `--runTreeShrink` to `--treeShrink`, `--consCatCutoff` to `--rateCutoff`, `--probaArg` to `--proba`, `--uncertaintyCoef` to `--phyloCoef`, `--internal` to `--internalBranches`, `--terminal` to `--terminalBranches`, and `--branchSpectra` to `--perBranch`.

### Fixed

- A nucleotide run no longer stops at the start because the taxonomy output is empty. `parsed_taxonomy.txt` is published only for protein input.
- The nucleotide outgroup id is matched case-sensitively, both when the header is checked and when that record is renamed. The id or species name has to be passed as it appears in the FASTA header.
- Deduplication moves the outgroup record to the front and then drops later exact copies of that sequence, so an identical ingroup sequence does not replace it. A reverse complement is kept.
- A gene that does not produce a synonymous spectrum is left out of `spectra_total.tsv` instead of adding an empty block.
- `--plot` defaults to true, the same value used in a normal run.
- `--minMuts` must be a whole number of at least 1, so `0` is rejected at startup rather than letting a header-only observation file through.
- Writing an ancestral-state table without its two input paths stops with a usage message.
- A TreeShrink crash no longer drops the gene. The run continues with the unpruned tree.

## [1.2.0] - 2026-09-07

### Added

- `--minMuts` (default `5`) sets how many reconstructed mutations a gene needs before a spectrum is derived. The old check compared the observation file’s size in bytes, so real spectra with few mutations were dropped and a large file could pass with almost no mutations.

### Fixed

- Nucleotide outgroups are taken from FASTA headers and renamed from that record, including when OUTGRP is not the last sequence and when encoded ids such as `seq_1` and `seq_10` would previously be rewritten together.
- A missing outgroup is reported and the protein path continues.
- Coding alignments drop gappy columns in whole codons so the reading frame stays intact.
- A short outgroup is kept through the length filter, and sequences are uppercased before encoding.
- A nucleotide `--input` that matches no files fails immediately.
- The high-N warning reads the N percentage.
- Spectrum calculation receives `--proba` when mutation extraction used probabilities.
- A BLAST database is accepted when `.ndb`, `.nin`, or `.nal` is present.

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
