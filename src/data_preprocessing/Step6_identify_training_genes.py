import pandas as pd
import glob
import sys
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
if CONFIG_DIR not in sys.path:
	sys.path.insert(0, CONFIG_DIR)
import config

# Define absolute file paths
# By default, read the least-likely genes from the intermediate file
# produced by Step4 under any '*/data_preprocessing/input/sampling/' folder.
# Fallback to config.least_likely_gene_path if not found.
search_pattern = os.path.join(ROOT_DIR, '**', 'data_preprocessing', 'input', 'sampling', 'least_likely_intermediate.tsv')
matches = glob.glob(search_pattern, recursive=True)
if matches:
	least_likely_file = matches[0]
	print(f"Using intermediate least-likely genes file: {least_likely_file}")
else:
	least_likely_file = config.least_likely_gene_path
	print(f"No intermediate least-likely file found; falling back to config path: {least_likely_file}")
#probable_genes_file = config.probable_gene_path
most_likely_genes_file = config.most_likely_gene_path
output_file = config.training_genes
ot_drugs = config.ot_phenotype_drugs

annoted_data_path = config.variant_output_directory

gwas_files = glob.glob(os.path.join(annoted_data_path, "Annotated_GWAS_*.csv"))
gwas_gene_dfs = [pd.read_csv(f)[['Gene']].drop_duplicates() for f in gwas_files]
gwas_genes_df = pd.concat(gwas_gene_dfs).drop_duplicates()

# --- Load full annotated GWAS files and compute genes that have at least one P<0.01 in EACH annotated file ---
# Requirement: a probable gene is retained only if it has at least one significant SNP (P < 0.01)
# in every Annotated_GWAS_*.csv file (intersection across per-file gene sets).
annot_pattern = os.path.join(config.variant_output_directory, "Annotated_GWAS_*.csv")
annot_files = sorted(glob.glob(annot_pattern))
if not annot_files:
	# fallback to single annotated_gwas path if present in config
	if hasattr(config, 'annotated_gwas') and os.path.exists(config.annotated_gwas):
		annot_files = [config.annotated_gwas]

significant_genes_perfile = set()
if annot_files:
	perfile_sets = []
	for p in annot_files:
		try:
			dfp = pd.read_csv(p)
		except Exception:
			dfp = pd.read_csv(p, engine='python')
		dfp['P'] = pd.to_numeric(dfp.get('P', pd.Series([])), errors='coerce')
		if 'Gene' in dfp.columns:
			genes_in_file = set(dfp[dfp['P'].notna() & (dfp['P'] < 0.01)]['Gene'].unique())
			perfile_sets.append(genes_in_file)
	if perfile_sets:
		# intersection: genes present with at least one significant SNP in every file
		significant_genes_perfile = set.intersection(*perfile_sets)
		print(f"Identified {len(significant_genes_perfile)} genes with at least one P<0.01 SNP in EACH Annotated_GWAS file (per-file intersection)")
	else:
		print("Warning: no Gene column found in annotated GWAS files; skipping p-value based filtering for probable genes")
else:
	print("No Annotated_GWAS files found; skipping p-value based filtering for probable genes")


# Read OT drugs file and select only 'gene' column, then remove duplicates
ot_drugs_df = pd.read_csv(ot_drugs, sep='\t')
ot_drugs_df = ot_drugs_df[['symbol']].drop_duplicates()
ot_drugs_df = ot_drugs_df.rename(columns={'symbol': 'Gene'})

# Read most likely genes file and select only 'Gene' column, then remove duplicates
most_likely_genes_df = pd.read_csv(most_likely_genes_file, sep='\t')
most_likely_genes_df = most_likely_genes_df[['Gene']].drop_duplicates()

# Filter GWAS genes to those in OT drugs and not in most likely genes
probable_genes_df = gwas_genes_df[
	gwas_genes_df['Gene'].isin(ot_drugs_df['Gene'])
	& ~gwas_genes_df['Gene'].isin(most_likely_genes_df['Gene'])
]
# Apply the per-file intersection filter: keep genes with at least one P<0.01 in every annotated file
if significant_genes_perfile:
	before_count = len(probable_genes_df)
	probable_genes_df = probable_genes_df[probable_genes_df['Gene'].isin(significant_genes_perfile)].copy()
	after_count = len(probable_genes_df)
	print(f"Probable genes: {before_count} -> {after_count} after applying per-file P<0.01 intersection filter")
else:
	probable_genes_df = probable_genes_df.copy()
	print("No genes matched the per-file P<0.01 intersection criterion; probable genes left unchanged.")

# Defensive: ensure probable_genes_df is a DataFrame with a proper index so we can
# safely assign a new column even if it is empty.
if not isinstance(probable_genes_df, pd.DataFrame):
	probable_genes_df = pd.DataFrame(probable_genes_df)
probable_genes_df = probable_genes_df.reset_index(drop=True)

# Read other gene files
least_likely_df = pd.read_csv(least_likely_file, sep='\t').copy()
most_likely_genes_df = pd.read_csv(most_likely_genes_file, sep='\t').copy()

# Add 'label' column to each DataFrame (use .loc to avoid SettingWithCopyWarning)
least_likely_df.loc[:, 'label'] = 'least likely'
probable_genes_df.loc[:, 'label'] = 'probable'
most_likely_genes_df.loc[:, 'label'] = 'most likely'

# Concatenate DataFrames and remove duplicate genes
combined_df = pd.concat([least_likely_df, probable_genes_df, most_likely_genes_df], ignore_index=True)
if 'Gene' in combined_df.columns:
	# keep one row per Gene, preserve first occurrence (labels should be set already)
	combined_df = combined_df.drop_duplicates(subset=['Gene'])
else:
	combined_df = combined_df.drop_duplicates()

# Ensure we write a compact training table (Gene + label if present)
if 'label' in combined_df.columns:
	out_df = combined_df[['Gene', 'label']].copy()
else:
	out_df = combined_df[['Gene']].copy()

# Save to output file (training genes)
out_df.to_csv(output_file, sep='\t', index=False)

# Also save the probable genes list to the configured path (if defined) to keep
# files in-sync. Write a deduplicated Gene-only table if possible.
try:
	probable_out = config.probable_gene_path
	if 'Gene' in probable_genes_df.columns:
		probable_genes_df[['Gene']].drop_duplicates().to_csv(probable_out, sep='\t', index=False)
	else:
		probable_genes_df.drop_duplicates().to_csv(probable_out, sep='\t', index=False)
except Exception:
	# if config does not define probable_gene_path or write fails, ignore
	pass
