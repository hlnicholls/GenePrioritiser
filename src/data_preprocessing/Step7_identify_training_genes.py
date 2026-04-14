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

# Define absolute file paths.
# Prefer pipeline-generated intermediate least-likely files in this repository,
# and avoid accidentally picking bundled example files.
preferred_candidates = [
	os.path.join(ROOT_DIR, 'input', 'training', 'sampling', 'least_likely_intermediate.tsv'),
	os.path.join(ROOT_DIR, 'results', 'data_preprocessing', 'input', 'sampling', 'least_likely_intermediate.tsv')
]

least_likely_file = None
for candidate in preferred_candidates:
	if os.path.exists(candidate):
		least_likely_file = candidate
		break

if least_likely_file is None:
	search_pattern = os.path.join(ROOT_DIR, '**', 'least_likely_intermediate.tsv')
	matches = glob.glob(search_pattern, recursive=True)
	# Exclude example/reference datasets.
	matches = [
		m for m in matches
		if '/example/' not in m and '/GenePrioritiser-db/' not in m
	]
	if matches:
		least_likely_file = matches[0]

if least_likely_file is None:
	least_likely_file = config.least_likely_gene_path
	print(f"No intermediate least-likely file found; falling back to config path: {least_likely_file}")
else:
	print(f"Using intermediate least-likely genes file: {least_likely_file}")
#probable_genes_file = config.probable_gene_path
most_likely_genes_file = config.most_likely_gene_path
output_file = config.training_genes
ot_drugs = config.ot_phenotype_drugs

annoted_data_path = config.variant_output_directory

# Read probable genes produced earlier by Step4. The selection of probable genes
# (OT drug overlap + per-file P<0.01 intersection) is performed in
# `Step4_probable_gene_selection.py`. Here we simply load the produced file.
try:
	probable_genes_df = pd.read_csv(config.probable_gene_path, sep='\t')
	# ensure a 'Gene' column
	if 'Gene' not in probable_genes_df.columns and probable_genes_df.shape[1] >= 1:
		probable_genes_df.columns = ['Gene'] + list(probable_genes_df.columns[1:])
	probable_genes_df = probable_genes_df[['Gene']].drop_duplicates().reset_index(drop=True)
	print(f"Loaded {len(probable_genes_df)} probable genes from: {config.probable_gene_path}")
except Exception as e:
	print(f"Failed to read probable genes from {config.probable_gene_path}: {e}")
	probable_genes_df = pd.DataFrame(columns=['Gene'])

# Defensive: ensure probable_genes_df is a DataFrame with a proper index so we can
# safely assign a new column even if it is empty.
if not isinstance(probable_genes_df, pd.DataFrame):
	probable_genes_df = pd.DataFrame(probable_genes_df)
probable_genes_df = probable_genes_df.reset_index(drop=True)

# Read other gene files
least_likely_df = pd.read_csv(least_likely_file, sep='\t').copy()
most_likely_genes_df = pd.read_csv(most_likely_genes_file, sep='\t').copy()
try:
	if 'Gene' in most_likely_genes_df.columns:
		ml_count = int(most_likely_genes_df['Gene'].drop_duplicates().shape[0])
	else:
		ml_count = int(most_likely_genes_df.shape[0])
	print(f"Loaded {ml_count} most likely genes from: {most_likely_genes_file}")
except Exception:
	print(f"Loaded most likely genes from: {most_likely_genes_file}")

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
out_dir = os.path.dirname(output_file)
if out_dir:
	os.makedirs(out_dir, exist_ok=True)
out_df.to_csv(output_file, sep='\t', index=False)

# Also save the probable genes list to the configured path (if defined) to keep
# files in-sync. Write a deduplicated Gene-only table if possible.
try:
	probable_out = config.probable_gene_path
	probable_dir = os.path.dirname(probable_out)
	if probable_dir:
		os.makedirs(probable_dir, exist_ok=True)
	if 'Gene' in probable_genes_df.columns:
		probable_genes_df[['Gene']].drop_duplicates().to_csv(probable_out, sep='\t', index=False)
	else:
		probable_genes_df.drop_duplicates().to_csv(probable_out, sep='\t', index=False)
except Exception:
	# if config does not define probable_gene_path or write fails, ignore
	pass
