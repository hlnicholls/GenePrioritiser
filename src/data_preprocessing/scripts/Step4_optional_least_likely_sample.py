import pandas as pd
import random
import sys
import os
CURRENT_DIR = os.getcwd()
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR, "..", "..", ".."))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
sys.path.append(CONFIG_DIR)
import config

# File paths
least_likely_path = config.least_likely_gene_path
most_likely_path = config.most_likely_gene_path
probable_path = config.probable_gene_path

# Read the files
least_likely_genes = pd.read_csv(least_likely_path, sep='\t')
most_likely_genes = pd.read_csv(most_likely_path, sep='\t')
probable_genes = pd.read_csv(probable_path, sep='\t')

# Calculate lengths
len_least_likely = len(least_likely_genes)
len_most_likely = len(most_likely_genes)
len_probable = len(probable_genes)
combined_length = len_most_likely + len_probable

# Check condition and sample if necessary
if len_least_likely > 2 * combined_length:
    sample_size = combined_length
    sampled_genes = least_likely_genes.sample(n=sample_size, random_state=1)
    sampled_genes.to_csv(least_likely_path, sep='\t', index=False)