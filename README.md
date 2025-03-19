# One_Health_OneTube

This repository contains scripts and concatenated metadata files used in the paper:

_Anker, K.M., Krog, J.S., Ciucani, M.M. and Trebbien, R._**"One Health in one tube: Optimising the workflow for whole genome sequencing of Influenza A viruses of human, swine, and avian origin"**, 2025


## Data
The repository includes the following concatenated metadata files used to compare sequencing outputs:

- `OT1_OT2_OT5_gel_scores.txt`
- `READ_COUNTS_all.txt`
- `coverage_all.txt.zip`
- `mapping_stats_all.txt`
- `passing_segments_all.txt`
- `sample_cts.txt`


## Scripts

#### IRMA Processing & Metadata Extraction:
Scripts to run [IRMA v1.0.2](https://doi.org/10.1186/s12864-016-3030-6) on raw sequencing data, extract metadata ([IRMA_postprocessing repository](https://github.com/marciux18/IRMA_postprocessing/tree/main)) and concatenate results from each sample.

#### Figure Generation
R scripts for plotting figures used in the paper.
