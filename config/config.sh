#!/usr/bin/env bash
# Shell configuration counterpart to config/config.py
# Define PROJECT_DIR and GWAS_PATHS here so shell scripts can source it.

PROJECT_DIR="${projectDir:-.}"

# List of GWAS files to process (example defaults). Users should edit this
# file to point to their GWAS summary statistics (GRCh37/hg19, gzipped).
GWAS_PATHS=(
  "$PROJECT_DIR/example/data_preprocessing/input/Evangelou_30224653_DBP.txt.gz"
  "$PROJECT_DIR/example/data_preprocessing/input/Evangelou_30224653_SBP.txt.gz"
  "$PROJECT_DIR/example/data_preprocessing/input/Evangelou_30224653_PP.txt.gz"
)

## Optionally override GTF path here
#GTF_PATH="$PROJECT_DIR/utils/gencode.v19.annotation.gtf.gz"

export PROJECT_DIR GWAS_PATHS
