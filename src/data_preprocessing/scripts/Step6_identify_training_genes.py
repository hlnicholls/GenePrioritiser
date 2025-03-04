import pandas as pd
import glob
import sys
import os
CURRENT_DIR = os.getcwd()
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
sys.path.append(CONFIG_DIR)
import config

# Define absolute file paths
least_likely_file = config.least_likely_gene_path
#probable_genes_file = config.probable_gene_path
most_likely_genes_file = config.most_likely_gene_path
output_file = config.training_genes
ot_drugs = config.ot_phenotype_drugs

annoted_data_path = config.variant_output_directory

gwas_files = glob.glob(os.path.join(annoted_data_path, "Annotated_GWAS_*.csv"))
gwas_gene_dfs = [pd.read_csv(f)[['Gene']].drop_duplicates() for f in gwas_files]
gwas_genes_df = pd.concat(gwas_gene_dfs).drop_duplicates()


# Read OT drugs file and select only 'gene' column, then remove duplicates
ot_drugs_df = pd.read_csv(ot_drugs, sep='\t')
ot_drugs_df = ot_drugs_df[['symbol']].drop_duplicates()
ot_drugs_df = ot_drugs_df.rename(columns={'symbol': 'Gene'})

# Read most likely genes file and select only 'Gene' column, then remove duplicates
most_likely_genes_df = pd.read_csv(most_likely_genes_file, sep='\t')
most_likely_genes_df = most_likely_genes_df[['Gene']].drop_duplicates()

# Filter GWAS genes to those in OT drugs and not in most likely genes
probable_genes_df = gwas_genes_df[gwas_genes_df['Gene'].isin(ot_drugs_df['Gene']) & ~gwas_genes_df['Gene'].isin(most_likely_genes_df['Gene'])]

# Read other gene files
least_likely_df = pd.read_csv(least_likely_file, sep='\t')
most_likely_genes_df = pd.read_csv(most_likely_genes_file, sep='\t')

# Add 'label' column to each DataFrame
least_likely_df['label'] = 'least likely'
probable_genes_df['label'] = 'probable'
most_likely_genes_df['label'] = 'most likely'

# Concatenate DataFrames
combined_df = pd.concat([least_likely_df, probable_genes_df, most_likely_genes_df], ignore_index=True)

# Save to output file
combined_df.to_csv(output_file, sep='\t', index=False)
