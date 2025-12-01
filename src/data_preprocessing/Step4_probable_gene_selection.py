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
		# find Gene column case-insensitively
		gene_col = None
		for c in dfp.columns:
			if str(c).strip().lower() == 'gene':
				gene_col = c
				break
		if gene_col is not None:
			g = pd.DataFrame(dfp[gene_col].drop_duplicates())
			g.columns = ['Gene']
			gwas_gene_frames.append(g)
	if gwas_gene_frames:
		gwas_genes_df = pd.concat(gwas_gene_frames).drop_duplicates().reset_index(drop=True)

# Determine intersection of genes that have at least one P < 0.01 in every file
significant_genes_perfile = set()
# p-value threshold for probable gene selection (can be overridden in config)
pv_thresh = getattr(config, 'probable_gene_pvalue_threshold', 0.01)
# Control whether genes must have a significant SNP in every file (intersection)
# or in any file (union). Default: union (at least one file) which matches
# the typical "gene needs at least one SNP" interpretation.
intersection_mode = getattr(config, 'probable_gene_intersection', False)
if annot_files:
	perfile_sets = []
	for p in annot_files:
		try:
			dfp = pd.read_csv(p)
		except Exception:
			dfp = pd.read_csv(p, engine='python')

		# detect gene column name
		gene_col = None
		for c in dfp.columns:
			if str(c).strip().lower() == 'gene':
				gene_col = c
				break

		# detect p-value column name (common variants)
		p_col = None
		for candidate in ['p', 'p_value', 'pvalue', 'p_val', 'pval', 'p.value', 'P']:
			for c in dfp.columns:
				if str(c).strip().lower() == candidate.lower():
					p_col = c
					break
			if p_col is not None:
				break

		if p_col is None:
			print(f"Warning: no p-value column found in {p}; this file will be ignored for p-value filtering")
			continue

		# coerce p-values
		dfp[p_col] = pd.to_numeric(dfp[p_col], errors='coerce')

		if gene_col is not None:
			genes_in_file = set(dfp[dfp[p_col].notna() & (dfp[p_col] < float(pv_thresh))][gene_col].dropna().astype(str).unique())
			perfile_sets.append(genes_in_file)
		else:
			print(f"Warning: 'Gene' column not found in {p}; this file will be ignored for p-value intersection")
	if perfile_sets:
		if intersection_mode:
			significant_genes_perfile = set.intersection(*perfile_sets)
			print(f"Identified {len(significant_genes_perfile)} genes with at least one P<{pv_thresh} SNP in EACH Annotated_GWAS file (per-file INTERSECTION)")
		else:
			significant_genes_perfile = set().union(*perfile_sets)
			print(f"Identified {len(significant_genes_perfile)} genes with at least one P<{pv_thresh} SNP in ANY Annotated_GWAS file (per-file UNION)")
	else:
		print("Warning: no valid Annotated_GWAS files with a 'Gene' column found; skipping p-value based filtering for probable genes")
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

# Apply per-file p-value intersection filter if annot files were present (even if intersection is empty)
if annot_files and 'perfile_sets' in locals():
	before_count = len(probable_genes_df)
	probable_genes_df = probable_genes_df[probable_genes_df['Gene'].isin(significant_genes_perfile)].copy()
	after_count = len(probable_genes_df)
	print(f"Probable genes: {before_count} -> {after_count} after applying per-file P<{pv_thresh} intersection filter")

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

