import pandas as pd
from itertools import chain
import sys
import os
CURRENT_DIR = os.getcwd()
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
sys.path.append(CONFIG_DIR)
import config


protein_alias = pd.read_csv(
    os.path.join(config.database_string_path , '9606.protein.aliases.v12.0.txt.gz'),
    sep='\t', compression='gzip', usecols=['#string_protein_id', 'alias']
)

protein_info = pd.read_csv(
    os.path.join(config.database_string_path , '9606.protein.info.v12.0.txt.gz'),
    sep='\t', compression='gzip', usecols=['#string_protein_id', 'preferred_name']
)

protein_links = pd.read_csv(
    os.path.join(config.database_string_path , '9606.protein.links.full.v12.0.txt.gz'),
    sep=' ', compression='gzip', usecols=['protein1', 'protein2', 'coexpression', 'experiments', 'database']
)

# Rename columns for consistency
protein_alias.columns = ['protein_id', 'Gene']
protein_info.columns = ['protein_id', 'Gene']

# Combine alias and info data, prioritizing preferred names from protein_info
combined_alias = pd.concat([protein_info, protein_alias]).drop_duplicates(subset='protein_id', keep='first')
alias_to_gene = combined_alias.set_index('protein_id')['Gene'].to_dict()

# Map protein IDs to gene names in protein_links
protein_links['Gene1'] = protein_links['protein1'].map(alias_to_gene).fillna('Unmapped')
protein_links['Gene2'] = protein_links['protein2'].map(alias_to_gene).fillna('Unmapped')

# Debug output for unmapped genes
unmapped_genes = protein_links[(protein_links['Gene1'] == 'Unmapped') | (protein_links['Gene2'] == 'Unmapped')]
print("Number of unmapped entries:", len(unmapped_genes))
print("Sample unmapped entries:")
print(unmapped_genes.head(10))

# Remove unmapped rows
protein_links = protein_links[(protein_links['Gene1'] != 'Unmapped') & (protein_links['Gene2'] != 'Unmapped')]

# Create sets for interactions
interaction_dict = {}
for _, row in protein_links.iterrows():
    interaction_dict.setdefault(row['Gene1'], set()).add(row['Gene2'])
    interaction_dict.setdefault(row['Gene2'], set()).add(row['Gene1'])

# Load known "most likely" genes
most_likely_genes = pd.read_csv(config.most_likely_gene_path, sep='\t')
most_likely_genes = set(most_likely_genes['Gene'].unique())

# Identify direct and indirect interactors for "most likely" genes
direct_and_indirect_interactors = set()
total_direct_interactors = set()
total_indirect_interactors = set()

for gene in most_likely_genes:
    if gene in interaction_dict:
        # Direct interactions
        direct_interactors = interaction_dict[gene]
        total_direct_interactors.update(direct_interactors)
        direct_and_indirect_interactors.update(direct_interactors)
        
        # Indirect interactions (one level away)
        for partner_gene in direct_interactors:
            if partner_gene in interaction_dict:
                indirect_interactors = interaction_dict[partner_gene]
                total_indirect_interactors.update(indirect_interactors)
                direct_and_indirect_interactors.update(indirect_interactors)

print("Number of direct interactors found with most likely disease genes:", len(total_direct_interactors))
print("Number of indirect interactors found with most likely disease genes:", len(total_indirect_interactors))

# Load and filter GWAS data based on P-value criteria
gwas_df = pd.read_csv(config.annotated_gwas)
filtered_gwas = gwas_df[gwas_df['P'] > 0.05]
valid_genes = filtered_gwas.groupby('Gene').filter(lambda x: all(x['P'] > 0.05))['Gene'].unique()

# Load gene types and filter for protein-coding genes
gene_types_df = pd.read_csv(config.gene_types, sep='\t', header=None)
gene_types_df.columns = ['chromosome', 'start', 'end', 'Gene', 'gene_type', 'source']
protein_coding_genes = set(gene_types_df[gene_types_df['gene_type'] == 'protein_coding']['Gene'].unique())

# Further filter valid genes to keep only protein-coding genes
valid_genes = [gene for gene in valid_genes if gene in protein_coding_genes]

# Filter valid genes not interacting with "most likely" or their secondary interactors
least_likely_genes = [gene for gene in valid_genes if gene not in direct_and_indirect_interactors]

# Save least likely genes
least_likely_genes_df = pd.DataFrame({'Gene': least_likely_genes})
least_likely_genes_df.to_csv(config.least_likely_gene_path, sep='\t', index=False)

if hasattr(config, "loci_ld_file") and os.path.exists(config.loci_ld_file):
    loci_ld_df = pd.read_csv(config.loci_ld_file, sep='\t', header=None)
    loci_ld_genes = set(loci_ld_df[0].unique())
    initial_gene_count = len(least_likely_genes)
    least_likely_genes = [gene for gene in least_likely_genes if gene not in loci_ld_genes]
    print(f"Filtered out {initial_gene_count - len(least_likely_genes)} genes based on Loci_LD.txt")

# Save least likely genes
least_likely_genes_df = pd.DataFrame({'Gene': least_likely_genes})
least_likely_genes_df.to_csv(config.least_likely_gene_path, sep='\t', index=False)


print("Least likely genes identified and saved.")
print("Number of least likely genes:", len(least_likely_genes))
print("Number of most likely genes:", len(most_likely_genes))
