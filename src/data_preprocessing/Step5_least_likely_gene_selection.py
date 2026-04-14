"""
Step3_least_likely_gene_selection.py

Purpose:
    Identify the "least-likely" genes for downstream negative-class construction.

What this script does (high level):
    - Load STRINGdb interaction files (aliases, preferred names, and full links).
    - Map protein IDs to gene symbols and build an undirected interaction adjacency map.
    - Load the "most likely" and optional "probable" gene lists (these are treated
        as positive seeds).
    - Perform a breadth-first search (BFS) from the seed genes to collect all genes
        within N hops. This is controlled by config.least_likely_ppi_mode:
        "direct" (1 hop) or "direct_and_secondary" (2 hops).
        Backward-compatible fallback: config.interactor_max_depth.
        These seed genes and their N-hop interactors are excluded from the negative set.
    - Load all Annotated_GWAS_*.csv files from the configured variant output
        directory, coerce P to numeric, and keep only rows with non-significant
        p-values (P > 0.05).
    - A gene is considered a candidate "least-likely" gene if ALL its annotated
        SNP rows are non-significant (P > 0.05) and it is protein coding and not
        within the seed/interactor neighborhood.
    - Optionally apply an extra filter file (config.least_likely_extra_filter).
    - Write the final least-likely gene list to config.least_likely_gene_path.

Inputs (from config):
    - config.database_string_path: location of STRINGdb files
    - config.most_likely_gene_path: path to most-likely genes (tsv, column 'Gene')
    - config.probable_gene_path (optional): path to probable genes (tsv)
    - config.variant_output_directory: directory containing Annotated_GWAS_*.csv
    - config.gene_types: local gene table used to restrict to protein_coding genes
    - config.least_likely_extra_filter (optional): extra gene exclusion list
    - config.least_likely_gene_path: where the resulting least-likely TSV is saved

Notes and assumptions:
    - Annotated GWAS files are expected to include a 'Gene' column and a 'P' column.
    - STRINGdb files are read from compressed TSVs; unmapped protein entries are
        removed before building interactions.
    - The interactor collection is BFS-based and includes genes reachable within
        max_depth hops (default 2 if no config is set). The seeds themselves are
        also excluded.
    - The output is a single-column TSV with header 'Gene'. This file is used by
        downstream sampling/ML steps.

This file is part of the data_preprocessing pipeline. Keep the behavior stable
because later steps expect the least-likely file to exist at the configured path.
"""

import pandas as pd
from itertools import chain
import sys
import os
import glob

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
if CONFIG_DIR not in sys.path:
    sys.path.insert(0, CONFIG_DIR)
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

# Optionally load "probable" genes (they should be excluded similarly)
probable_genes = set()
if hasattr(config, 'probable_gene_path') and os.path.exists(config.probable_gene_path):
    try:
        pg = pd.read_csv(config.probable_gene_path, sep='\t')
        probable_genes = set(pg['Gene'].unique())
    except Exception:
        print(f"Warning: could not read probable genes from {config.probable_gene_path}")

# Collect interactors up to a configurable depth (N hops) using BFS
from collections import deque

def collect_interactors(seed_genes, interaction_dict, max_depth=2):
    """Return set of genes reachable from any seed within max_depth hops.
    max_depth=1 -> direct interactors; max_depth=2 -> direct + one-level indirect, etc.
    """
    seen = set()
    result = set()
    q = deque()
    # initialize queue: (gene, depth)
    for g in seed_genes:
        seen.add(g)
        q.append((g, 0))

    while q:
        gene, depth = q.popleft()
        if depth >= max_depth:
            continue
        neighbors = interaction_dict.get(gene, ())
        for nb in neighbors:
            if nb not in seen:
                seen.add(nb)
                q.append((nb, depth + 1))
                # neighbor is within 1..max_depth of a seed
                result.add(nb)
    return result


MAX_INTERACTOR_DEPTH = getattr(config, 'interactor_max_depth', 2)

# Build seed set from most-likely and probable genes
seed_genes = set(most_likely_genes)
if probable_genes:
    seed_genes = seed_genes.union(probable_genes)

# Resolve PPI filtering mode.
# Supported modes (config.least_likely_ppi_mode):
#   - "direct"                -> exclude only direct interactors (depth=1)
#   - "direct_and_secondary"  -> exclude direct + second-degree interactors (depth=2)
# Backward compatibility:
#   - If least_likely_ppi_mode is not set, fall back to config.interactor_max_depth (if present),
#     else default to direct_and_secondary (depth=2).
ppi_mode = getattr(config, 'least_likely_ppi_mode', None)
if ppi_mode is not None:
    ppi_mode = str(ppi_mode).strip().lower()

if ppi_mode in {'direct', 'direct_and_secondary'}:
    MAX_INTERACTOR_DEPTH = 1 if ppi_mode == 'direct' else 2
    print(f"Using PPI filter mode: {ppi_mode} (max_depth={MAX_INTERACTOR_DEPTH})")
else:
    # legacy behavior: explicit interactor_max_depth still works
    MAX_INTERACTOR_DEPTH = getattr(config, 'interactor_max_depth', 2)
    if ppi_mode is not None:
        print(
            f"Warning: unsupported least_likely_ppi_mode='{ppi_mode}'. "
            f"Falling back to interactor_max_depth={MAX_INTERACTOR_DEPTH}."
        )
    else:
        print(
            "No least_likely_ppi_mode set; using interactor_max_depth="
            f"{MAX_INTERACTOR_DEPTH} (legacy fallback)."
        )

# Collect interactors up to MAX_INTERACTOR_DEPTH
direct_and_indirect_interactors = collect_interactors(
    seed_genes,
    interaction_dict,
    max_depth=MAX_INTERACTOR_DEPTH,
)

# For reporting, also collect stats for direct only and one-level for diagnostics
direct_only = collect_interactors(seed_genes, interaction_dict, max_depth=1)
one_level = collect_interactors(seed_genes, interaction_dict, max_depth=2)

print(f"Interactors collected with max_depth={MAX_INTERACTOR_DEPTH}: {len(direct_and_indirect_interactors)} genes")
print(f"Direct only interactors (depth=1): {len(direct_only)}")
print(f"One-level interactors (depth=2): {len(one_level)}")

# Also exclude the seed genes themselves from being labelled as least-likely
direct_and_indirect_interactors = direct_and_indirect_interactors.union(seed_genes)

pattern = os.path.join(config.variant_output_directory, "Annotated_GWAS_*.csv")
annot_files = sorted(glob.glob(pattern))
if not annot_files:
    # fallback: if a single annotated_gwas path is configured use that
    if hasattr(config, 'annotated_gwas') and os.path.exists(config.annotated_gwas):
        annot_files = [config.annotated_gwas]
    else:
        raise FileNotFoundError(f"No annotated GWAS files found with pattern: {pattern}")

print(f"Found {len(annot_files)} annotated GWAS file(s) to load:")
for p in annot_files:
    print(' -', p)

df_list = []
for p in annot_files:
    try:
        df_list.append(pd.read_csv(p))
    except Exception:
        # try with automatic engine detection (some files might be gz or different separators)
        df_list.append(pd.read_csv(p, engine='python'))

gwas_df = pd.concat(df_list, ignore_index=True)

# detect gene column (case-insensitive)
gene_col = None
for c in gwas_df.columns:
    if str(c).strip().lower() == 'gene':
        gene_col = c
        break
if gene_col is None:
    raise KeyError('Gene column not found in annotated GWAS files')

# detect p-value column (case-insensitive)
p_col = None
for candidate in ['p', 'p_value', 'pvalue', 'p_val', 'pval', 'p.value', 'P']:
    for c in gwas_df.columns:
        if str(c).strip().lower() == candidate.lower():
            p_col = c
            break
    if p_col is not None:
        break
if p_col is None:
    raise KeyError('No p-value column found in annotated GWAS files')

# p-value threshold for least-likely selection (can be overridden in config)
pv_thresh = getattr(config, 'least_likely_pvalue_threshold', 0.05)

print(f"Using gene column: '{gene_col}', p-value column: '{p_col}', threshold: {pv_thresh}")

# Ensure p-values are numeric
gwas_df[p_col] = pd.to_numeric(gwas_df[p_col], errors='coerce')

# Select genes where ALL rows for that gene are non-significant (p > pv_thresh)
# and p-values are present (not NaN)
grouped = gwas_df.groupby(gene_col)
valid_genes = []
for g, grp in grouped:
    grp_p = grp[p_col]
    if grp_p.notna().all() and (grp_p > pv_thresh).all():
        valid_genes.append(str(g))
valid_genes = pd.unique(valid_genes)

# Load gene types and filter for protein-coding genes
gene_types_df = pd.read_csv(config.gene_types, sep='\t', header=None)
gene_types_df.columns = ['chromosome', 'start', 'end', 'Gene', 'gene_type', 'source']
protein_coding_genes = set(gene_types_df[gene_types_df['gene_type'] == 'protein_coding']['Gene'].unique())

# Further filter valid genes to keep only protein-coding genes
valid_genes = [gene for gene in valid_genes if gene in protein_coding_genes]

# Filter valid genes not interacting with "most likely" or their secondary interactors
least_likely_genes = [gene for gene in valid_genes if gene not in direct_and_indirect_interactors]

# Save least likely genes to an intermediate sampling file (do not overwrite
# the final least_likely_genes path here). Place this intermediate file next
# to the configured final least_likely path (under a 'sampling' subdir) so
# it follows the user's configured PROJECT_DIR (e.g., 'results').
out_dir = os.path.join(os.path.dirname(config.least_likely_gene_path), 'sampling')
os.makedirs(out_dir, exist_ok=True)
intermediate_least_path = os.path.join(out_dir, 'least_likely_intermediate.tsv')
least_likely_genes_df = pd.DataFrame({'Gene': least_likely_genes})
least_likely_genes_df.to_csv(intermediate_least_path, sep='\t', index=False)
print(f'Wrote intermediate least-likely genes to {intermediate_least_path} (will not overwrite {config.least_likely_gene_path})')

if hasattr(config, "least_likely_extra_filter") and os.path.exists(config.least_likely_extra_filter):
    print("Applying additional filtering based on:", config.least_likely_extra_filter)
    extra_df = pd.read_csv(config.least_likely_extra_filter, sep='\t', header=None)
    extra_genes = set(extra_df[0].unique())
    initial_gene_count = len(least_likely_genes)
    least_likely_genes = [gene for gene in least_likely_genes if gene not in extra_genes]
    print(f"Filtered out {initial_gene_count - len(least_likely_genes)} genes based on additional criteria.")

# Save least likely genes after optional additional filtering to the same
# intermediate file (do not overwrite final least_likely file here)
least_likely_genes_df = pd.DataFrame({'Gene': least_likely_genes})
least_likely_genes_df.to_csv(intermediate_least_path, sep='\t', index=False)
print(f'Updated intermediate least-likely genes at {intermediate_least_path}')


print("Least likely genes identified and saved.")
print("Number of least likely genes:", len(least_likely_genes))
print("Number of probable genes:", len(probable_genes))
print("Number of most likely genes:", len(most_likely_genes))
