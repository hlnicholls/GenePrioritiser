#!/usr/bin/env bash
# Shell configuration counterpart to config/config.py
# Define PROJECT_DIR and GWAS_PATHS here so shell scripts can source it.

PROJECT_DIR="${projectDir:-.}"

# List of GWAS files to process (example defaults). Users should edit this
# file to point to their GWAS summary statistics (GRCh37/hg19, gzipped).
GWAS_PATHS=(
  "$PROJECT_DIR/results/data_preprocessing/input/GCST90310294_SBP.tsv.gz"
  "$PROJECT_DIR/results/data_preprocessing/input/GCST90310295_DBP.tsv.gz"
  "$PROJECT_DIR/results/data_preprocessing/input/GCST90310296_PP.tsv.gz"
)

## Optionally override GTF path here
#GTF_PATH="$PROJECT_DIR/utils/gencode.v19.annotation.gtf.gz"

export PROJECT_DIR GWAS_PATHS
