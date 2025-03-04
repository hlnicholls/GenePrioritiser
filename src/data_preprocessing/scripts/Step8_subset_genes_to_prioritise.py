import pandas as pd
import sys
import os
CURRENT_DIR = os.getcwd()
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
sys.path.append(CONFIG_DIR)
import config

# Load datasets
all_genes_df = pd.read_csv(config.all_genes_all_features_unprocessed)
training_genes_df = pd.read_csv(config.training_genes_features)

# Identify genes to keep
# Genes in GWAS but not in training
genes_to_keep = set(all_genes_df["Gene"]) - set(training_genes_df["Gene"])

# Subset the all_genes_df to keep only the desired genes
subset_genes_df = all_genes_df[all_genes_df["Gene"].isin(genes_to_keep)]

# Save the result
subset_genes_df.to_csv(config.genes_to_prioritise, index=False)
