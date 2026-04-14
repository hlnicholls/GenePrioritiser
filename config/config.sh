#!/usr/bin/env bash
# Shell configuration counterpart to config/config.py
# Define PROJECT_DIR and GWAS_PATHS here so shell scripts can source it.

PROJECT_DIR="${projectDir:-.}"

# List of GWAS files to process (example defaults). Users should edit this
# file to point to their GWAS summary statistics (GRCh37/hg19, gzipped).
GWAS_PATHS=(
  "$PROJECT_DIR/input/full_gwas/SBP_GCST90310294.tsv.gz"
  "$PROJECT_DIR/input/full_gwas/DBP_GCST90310295.tsv.gz"
  "$PROJECT_DIR/input/full_gwas/PP_GCST90310296.tsv.gz"
)

## Optionally override GTF path here
#GTF_PATH="$PROJECT_DIR/utils/gencode.v19.annotation.gtf.gz"

export PROJECT_DIR GWAS_PATHS
