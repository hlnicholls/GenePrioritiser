import os
import glob
import sys
import pandas as pd
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
if CONFIG_DIR not in sys.path:
	sys.path.insert(0, CONFIG_DIR)
import config

# This script selects "probable" genes and writes them to the path defined in
# `config.probable_gene_path`. Criteria:
#  - Gene appears in at least one Annotated_GWAS_*.csv file
#  - Gene is present in OpenTargets drug list for the phenotype (config.ot_phenotype_drugs)
#  - Gene is NOT in the most likely gene set (config.most_likely_gene_path)
#  - Additionally, retain only genes that have at least one SNP with P < 0.01 in EVERY
#    Annotated_GWAS_*.csv file (intersection across files). If no annotated files are
#    present, the per-file P-value filter is skipped.

annot_pattern = os.path.join(config.variant_output_directory, "Annotated_GWAS_*.csv")
annot_files = sorted(glob.glob(annot_pattern))
if not annot_files:
	if hasattr(config, 'annotated_gwas') and os.path.exists(config.annotated_gwas):
		annot_files = [config.annotated_gwas]

# Collect all genes present in Annotated_GWAS files (if present)
gwas_gene_sets = []
gwas_genes_df = pd.DataFrame(columns=['Gene'])
if annot_files:
	gwas_gene_frames = []
	for p in annot_files:
		try:
			dfp = pd.read_csv(p)
		except Exception:
			dfp = pd.read_csv(p, engine='python')
		if 'Gene' in dfp.columns:
			g = pd.DataFrame(dfp['Gene'].drop_duplicates())
			g.columns = ['Gene']
			gwas_gene_frames.append(g)
	if gwas_gene_frames:
		gwas_genes_df = pd.concat(gwas_gene_frames).drop_duplicates().reset_index(drop=True)

# Determine intersection of genes that have at least one P < 0.01 in every file
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
		significant_genes_perfile = set.intersection(*perfile_sets)
		print(f"Identified {len(significant_genes_perfile)} genes with at least one P<0.01 SNP in EACH Annotated_GWAS file (per-file intersection)")
	else:
		print("Warning: no Gene column found in annotated GWAS files; skipping p-value based filtering for probable genes")
else:
	print("No Annotated_GWAS files found; probable gene selection will use OT-drugs overlap only")

# Read OpenTargets drugs for phenotype
try:
	ot_drugs_df = pd.read_csv(config.ot_phenotype_drugs, sep='\t')
	if 'symbol' in ot_drugs_df.columns:
		ot_drugs_df = ot_drugs_df[['symbol']].drop_duplicates().rename(columns={'symbol': 'Gene'})
	elif 'Gene' in ot_drugs_df.columns:
		ot_drugs_df = ot_drugs_df[['Gene']].drop_duplicates()
	else:
		ot_drugs_df = pd.DataFrame(columns=['Gene'])
except Exception as e:
	print(f"Failed to read OpenTargets drugs file {config.ot_phenotype_drugs}: {e}")
	ot_drugs_df = pd.DataFrame(columns=['Gene'])

# Read most likely genes
try:
	most_likely_genes_df = pd.read_csv(config.most_likely_gene_path, sep='\t')
	if 'Gene' in most_likely_genes_df.columns:
		most_likely_genes_df = most_likely_genes_df[['Gene']].drop_duplicates()
	else:
		most_likely_genes_df = pd.DataFrame(columns=['Gene'])
except Exception as e:
	print(f"Failed to read most likely genes file {config.most_likely_gene_path}: {e}")
	most_likely_genes_df = pd.DataFrame(columns=['Gene'])

# Start from GWAS genes if available, else try to use OT list as base
if not gwas_genes_df.empty:
	base_genes = gwas_genes_df['Gene']
else:
	base_genes = ot_drugs_df['Gene']

# Filter: in OT drug list and not in most likely genes
probable_genes_df = pd.DataFrame({'Gene': base_genes.unique()})
if not ot_drugs_df.empty:
	probable_genes_df = probable_genes_df[probable_genes_df['Gene'].isin(ot_drugs_df['Gene'])]
if not most_likely_genes_df.empty:
	probable_genes_df = probable_genes_df[~probable_genes_df['Gene'].isin(most_likely_genes_df['Gene'])]

# Apply per-file P<0.01 intersection filter if available
if significant_genes_perfile:
	before_count = len(probable_genes_df)
	probable_genes_df = probable_genes_df[probable_genes_df['Gene'].isin(significant_genes_perfile)].copy()
	after_count = len(probable_genes_df)
	print(f"Probable genes: {before_count} -> {after_count} after applying per-file P<0.01 intersection filter")

# Ensure Gene column exists and deduplicate
if 'Gene' in probable_genes_df.columns:
	probable_genes_df = probable_genes_df[['Gene']].drop_duplicates().reset_index(drop=True)
else:
	probable_genes_df = pd.DataFrame(columns=['Gene'])

# Write probable genes to configured path
try:
	out_path = config.probable_gene_path
	probable_genes_df.to_csv(out_path, sep='\t', index=False)
	print(f"Wrote {len(probable_genes_df)} probable genes to: {out_path}")
except Exception as e:
	print(f"Failed to write probable genes to {getattr(config, 'probable_gene_path', 'unknown')}: {e}")

