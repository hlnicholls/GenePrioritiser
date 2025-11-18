"""
Stratified undersampling of the 'least-likely' gene set using feature completeness.

This script is intended for the multiclass setup where you have three classes:
    - most-likely (positive)
    - probable  (positive-ish)
    - least-likely (negative)

New behaviour (run after the database merge step):
    - read merged feature table (all genes x features) produced by the merge step
        (config.all_genes_all_features_unprocessed)
    - compute per-gene completeness = fraction of non-missing feature values
    - bin genes by variant_count (derived from Annotated_GWAS files) using the
        combined positives distribution (as before)
    - within each bin, preferentially select least-likely genes that have the
        highest completeness (i.e. most complete features)
    - still preserve the overall balanced target size (match combined positives)

This makes the sampled negatives more usable for ML since they contain more
complete feature vectors.

Usage: run without args (reads paths from config). Optional args:
    --seed INT    reproducible seed (default 42)
    --n_bins INT  number of variant_count bins (default 3 - tertiles)
    --force       run sampling even if least <= combined

"""

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..'))
CONFIG_DIR = os.path.join(ROOT_DIR, 'config')
if CONFIG_DIR not in sys.path:
    sys.path.insert(0, CONFIG_DIR)
import config


def _detect_variant_col(df):
    candidates = ['SNP', 'snp', 'Marker', 'marker', 'MarkerName', 'rsid', 'rsID', 'ID', 'markername']
    for c in candidates:
        if c in df.columns:
            return c
    # fallback: first column that's not the gene-like columns
    for c in df.columns:
        if c.lower() not in ('gene', 'gene_symbol', 'gene symbol'):
            return c
    return df.columns[0]


def _detect_gene_col(df):
    for c in ('Gene', 'gene', 'Gene symbol', 'Gene_Symbol', 'gene_symbol'):
        if c in df.columns:
            return c
    # fallback to first column
    return df.columns[0]


def main(seed=42, n_bins=3, force=False):
    # Paths from config
    # By default, read the least-likely genes from the intermediate file
    # produced by the least-likely selection step. We search for the file
    # 'least_likely_intermediate.tsv' under any '*/data_preprocessing/input/sampling/'
    # directory starting from the repository root so the pipeline works even if
    # the parent (e.g. 'example' vs 'utils/example') changes.
    search_pattern = os.path.join(ROOT_DIR, '**', 'data_preprocessing', 'input', 'sampling', 'least_likely_intermediate.tsv')
    matches = glob.glob(search_pattern, recursive=True)
    if matches:
        least_path = matches[0]
        print(f"Using intermediate least-likely genes file: {least_path}")
    else:
        # fallback to config path if no intermediate file found
        least_path = config.least_likely_gene_path
        print(f"No intermediate least-likely file found; falling back to config path: {least_path}")
    most_path = config.most_likely_gene_path
    probable_path = config.probable_gene_path
    ann_dir = config.variant_output_directory

    annotated_glob = os.path.join(ann_dir, 'Annotated_GWAS_*.csv')
    ann_files = sorted(glob.glob(annotated_glob))
    if len(ann_files) == 0:
        raise SystemExit(f'No annotated GWAS files found at {annotated_glob}')

    # 1) compute variant counts per gene across annotated GWAS files
    parts = []
    for f in ann_files:
        df = pd.read_csv(f)
        var_col = _detect_variant_col(df)
        gene_col = _detect_gene_col(df)
        parts.append(df[[var_col, gene_col]].dropna().rename(columns={var_col: 'variant_id', gene_col: 'Gene'}))
    all_ann = pd.concat(parts, ignore_index=True)
    variant_counts = all_ann.groupby('Gene')['variant_id'].nunique().reset_index(name='variant_count')

    # 3) load class files
    least_df = pd.read_csv(least_path, sep='\t')
    most_df = pd.read_csv(most_path, sep='\t')
    probable_df = pd.read_csv(probable_path, sep='\t')

    # ensure a consistent gene column name
    for df in (least_df, most_df, probable_df):
        if 'Gene' not in df.columns and 'gene' in df.columns:
            df.rename(columns={'gene': 'Gene'}, inplace=True)
        if 'Gene' not in df.columns:
            df.rename(columns={df.columns[0]: 'Gene'}, inplace=True)

    # 3) merge variant_count into each set (missing -> 0)
    least_df = least_df.merge(variant_counts, on='Gene', how='left').fillna({'variant_count': 0})
    most_df = most_df.merge(variant_counts, on='Gene', how='left').fillna({'variant_count': 0})
    probable_df = probable_df.merge(variant_counts, on='Gene', how='left').fillna({'variant_count': 0})

    # 3b) read merged feature table to compute completeness for each gene.
    # This file is produced by the merge step and should exist when this
    # script is run after that step.
    completeness_col = 'completeness'
    try:
        features_df = pd.read_csv(config.all_genes_all_features_unprocessed)
        # ensure 'Gene' column exists
        if 'Gene' not in features_df.columns and 'gene' in features_df.columns:
            features_df.rename(columns={'gene': 'Gene'}, inplace=True)
        feature_cols = [c for c in features_df.columns if c != 'Gene']
        if len(feature_cols) == 0:
            # no features: set completeness to 0
            features_df[completeness_col] = 0.0
        else:
            features_df[completeness_col] = features_df[feature_cols].notnull().sum(axis=1) / float(len(feature_cols))
        features_df = features_df[['Gene', completeness_col]]
    except Exception as e:
        print('Warning: could not read features file to compute completeness:', e)
        # create an empty completeness DF so merge results in completeness=0
        features_df = pd.DataFrame(columns=['Gene', completeness_col])

    # merge completeness into least_df (and also into combined positives for diagnostics)
    least_df = least_df.merge(features_df, on='Gene', how='left').fillna({completeness_col: 0.0})
    most_df = most_df.merge(features_df, on='Gene', how='left').fillna({completeness_col: 0.0})
    probable_df = probable_df.merge(features_df, on='Gene', how='left').fillna({completeness_col: 0.0})

    # Use ONLY the most-likely gene set as the positive reference for sampling
    # (match previous behaviour: sample negatives to equal the most-likely size)
    combined_pos = most_df[['Gene', 'variant_count']].copy()
    combined_size = len(combined_pos)
    least_size = len(least_df)

    if combined_size == 0:
        raise SystemExit('Most-likely set size is zero; cannot sample to match.')

    # target_total is the size of the most-likely group
    target_total = combined_size

    if not force and least_size <= target_total:
        print(f'No sampling needed: least size {least_size} <= most-likely size {target_total} (use --force to override).')
        # Still update sampled/least files and regenerate the training files so
        # downstream steps see the (unchanged) least-likely set and joined data.
        sampled = least_df.copy()

        out_dir = os.path.join(ROOT_DIR, 'example', 'data_preprocessing', 'input', 'sampling')
        os.makedirs(out_dir, exist_ok=True)
        sampled_path = os.path.join(out_dir, f'least_likely_sampled_seed{seed}.tsv')
        audit_path = os.path.join(out_dir, f'least_likely_sampling_audit_seed{seed}.json')

        sampled[['Gene']].drop_duplicates().to_csv(sampled_path, sep='\t', index=False)
        sampled[['Gene']].drop_duplicates().to_csv(config.least_likely_gene_path, sep='\t', index=False)

        audit = {
            'seed': int(seed),
            'target_total': int(target_total),
            'original_least': int(least_size),
            'note': 'no sampling performed; least set unchanged'
        }
        try:
            with open(audit_path, 'w') as fh:
                json.dump(audit, fh, indent=2)
        except Exception:
            pass

        # write combined training_genes and joined feature files as if sampling had occurred
        try:
            most_df_small = most_df.copy()
            probable_df_small = probable_df.copy()

            def _ensure_gene_col(df):
                if 'Gene' not in df.columns and 'gene' in df.columns:
                    df = df.rename(columns={'gene': 'Gene'})
                if 'Gene' not in df.columns and df.shape[1] >= 1:
                    df = df.rename(columns={df.columns[0]: 'Gene'})
                return df[['Gene']].drop_duplicates()

            most_df_small = _ensure_gene_col(most_df_small)
            probable_df_small = _ensure_gene_col(probable_df_small)
            least_sampled = sampled[['Gene']].drop_duplicates()

            most_df_small['label'] = 'most likely'
            probable_df_small['label'] = 'probable'
            least_sampled['label'] = 'least likely'

            combined_df = pd.concat([least_sampled, probable_df_small, most_df_small], ignore_index=True)
            combined_df.drop_duplicates(subset=['Gene'], inplace=True)
            try:
                combined_df.to_csv(config.training_genes, sep='\t', index=False)
                print('Wrote updated training_genes to', config.training_genes)
            except Exception as e:
                print('Warning: could not write updated training_genes file:', e)

            # recreate joined training features and EDA input
            try:
                training_df = pd.read_csv(config.training_genes, sep='\t')
                ml_df = pd.read_csv(config.all_genes_all_features_unprocessed)
                result_df = pd.merge(training_df, ml_df, on='Gene', how='left')
                # drop duplicates by Gene so training table has one row per gene
                if 'Gene' in result_df.columns:
                    result_df.drop_duplicates(subset=['Gene'], inplace=True)
                else:
                    result_df.drop_duplicates(inplace=True)
                result_df.to_csv(config.training_genes_features, index=False)
                result_df.to_csv(config.training_data_all_features_eda, index=False)
                print('Wrote training features to:', config.training_genes_features)
                print('Wrote EDA training input to:', config.training_data_all_features_eda)
            except Exception as e:
                print('Warning: could not join training_genes with all features:', e)
        except Exception as e:
            print('Warning while writing unchanged sampled/training files:', e)

        return

    # 4) bin genes by variant_count using quantile-derived edges from combined positives
    # we derive edges from combined_pos to make the negative sampling distribution match positives
    quantiles = np.linspace(0, 1, n_bins + 1)
    edges = combined_pos['variant_count'].quantile(quantiles).values
    edges = np.unique(edges)
    if len(edges) <= 1:
        # degenerate: everyone has same count, create two bins [min, max+1]
        edges = np.array([combined_pos['variant_count'].min(), combined_pos['variant_count'].max() + 1])

    combined_pos['bin'] = pd.cut(combined_pos['variant_count'], bins=edges, include_lowest=True, labels=False)
    combined_pos['bin'] = combined_pos['bin'].fillna(0).astype(int)

    least_df['bin'] = pd.cut(least_df['variant_count'], bins=edges, include_lowest=True, labels=False)
    least_df['bin'] = least_df['bin'].fillna(0).astype(int)

    # 5) compute desired counts per bin proportional to combined_pos distribution
    pos_counts_by_bin = combined_pos['bin'].value_counts().sort_index()
    pos_props = (pos_counts_by_bin / pos_counts_by_bin.sum()).to_dict()
    # target_total already set to most-likely size above
    target_by_bin = {int(b): int(round(pos_props.get(b, 0) * target_total)) for b in pos_props.keys()}
    # correct rounding error
    cur_sum = sum(target_by_bin.values())
    if cur_sum != target_total and len(target_by_bin) > 0:
        diff = target_total - cur_sum
        largest_bin = max(target_by_bin, key=target_by_bin.get)
        target_by_bin[largest_bin] += diff

    # 6) sample from least per bin, prioritising genes with higher completeness
    np.random.seed(seed)
    sampled_parts = []
    audit = {
        'seed': int(seed),
        'target_total': int(target_total),
        'original_least': int(least_size),
        'per_bin': {}
    }
    for b, desired in target_by_bin.items():
        bin_df = least_df[least_df['bin'] == b].copy()
        available = len(bin_df)
        take = min(available, int(desired))
        audit['per_bin'][str(b)] = {'desired': int(desired), 'available': int(available), 'taken': int(take)}
        if take > 0 and available > 0:
            # sort by completeness desc, then variant_count desc as a tie-breaker
            # break ties randomly but deterministically using seed
            bin_df = bin_df.sample(frac=1, random_state=seed).sort_values([completeness_col, 'variant_count'], ascending=[False, False])
            sampled_parts.append(bin_df.head(take))

    sampled = pd.concat(sampled_parts, ignore_index=True) if len(sampled_parts) > 0 else pd.DataFrame(columns=least_df.columns)
    # top-up if some bins lacked enough genes
    if len(sampled) < target_total:
        remaining_needed = int(target_total - len(sampled))
        remaining_pool = least_df[~least_df['Gene'].isin(sampled['Gene'])].copy()
        # prefer most complete in remaining pool
        remaining_pool = remaining_pool.sample(frac=1, random_state=seed).sort_values([completeness_col, 'variant_count'], ascending=[False, False])
        take = min(len(remaining_pool), remaining_needed)
        if take > 0:
            sampled = pd.concat([sampled, remaining_pool.head(take)], ignore_index=True)
        audit['topped_up'] = int(take)

    # 7) create output directory under example/data_preprocessing/input/sampling
    out_dir = os.path.join(ROOT_DIR, 'example', 'data_preprocessing', 'input', 'sampling')
    os.makedirs(out_dir, exist_ok=True)
    # Create diagnostic plots: distribution of variant_count for combined positives and least
    try:
        plt.figure(figsize=(8, 4))
        sns.histplot(combined_pos['variant_count'], color='C0', label='combined positives', kde=False, stat='density', bins=50, alpha=0.6)
        sns.histplot(least_df['variant_count'], color='C1', label='least-likely', kde=False, stat='density', bins=50, alpha=0.6)
        plt.legend()
        plt.xlabel('variant_count')
        plt.title('Variant count distribution: combined positives vs least-likely')
        hist_path = os.path.join(out_dir, f'variant_count_hist_seed{seed}.png')
        plt.tight_layout()
        plt.savefig(hist_path, dpi=150)
        plt.close()

        plt.figure(figsize=(6, 4))
        sns.boxplot(data=[combined_pos['variant_count'], least_df['variant_count']], palette=['C0', 'C1'])
        plt.xticks([0, 1], ['combined_pos', 'least'])
        plt.ylabel('variant_count')
        plt.title('Variant count boxplot')
        box_path = os.path.join(out_dir, f'variant_count_box_seed{seed}.png')
        plt.tight_layout()
        plt.savefig(box_path, dpi=150)
        plt.close()
        print('Wrote diagnostic plots:', hist_path, box_path)
    except Exception as e:
        print('Warning: could not create diagnostic plots:', e)
    sampled_path = os.path.join(out_dir, f'least_likely_sampled_seed{seed}.tsv')
    audit_path = os.path.join(out_dir, f'least_likely_sampling_audit_seed{seed}.json')
    sampled[['Gene']].drop_duplicates().to_csv(sampled_path, sep='\t', index=False)
    sampled[['Gene']].drop_duplicates().to_csv(config.least_likely_gene_path, sep='\t', index=False)
    with open(audit_path, 'w') as fh:
        json.dump(audit, fh, indent=2)

    # Also update the combined training genes file used by the pipeline so that
    # downstream steps use the filtered least-likely set. This overwrites
    # config.training_genes with a combined table of Gene + label columns.
    try:
        # read most and probable lists if available
        most_df = pd.read_csv(most_path, sep='\t') if most_path and os.path.exists(most_path) else pd.DataFrame(columns=['Gene'])
    except Exception:
        most_df = pd.DataFrame(columns=['Gene'])
    try:
        probable_df = pd.read_csv(probable_path, sep='\t') if probable_path and os.path.exists(probable_path) else pd.DataFrame(columns=['Gene'])
    except Exception:
        probable_df = pd.DataFrame(columns=['Gene'])

    def _ensure_gene_col(df):
        if 'Gene' not in df.columns and 'gene' in df.columns:
            df = df.rename(columns={'gene': 'Gene'})
        if 'Gene' not in df.columns and df.shape[1] >= 1:
            df = df.rename(columns={df.columns[0]: 'Gene'})
        return df[['Gene']].drop_duplicates()

    most_df = _ensure_gene_col(most_df)
    probable_df = _ensure_gene_col(probable_df)
    least_sampled = sampled[['Gene']].drop_duplicates()

    most_df['label'] = 'most likely'
    probable_df['label'] = 'probable'
    least_sampled['label'] = 'least likely'

    combined_df = pd.concat([least_sampled, probable_df, most_df], ignore_index=True)
    combined_df.drop_duplicates(subset=['Gene'], inplace=True)
    try:
        combined_df.to_csv(config.training_genes, sep='\t', index=False)
        print('Wrote updated training_genes to', config.training_genes)
    except Exception as e:
        print('Warning: could not write updated training_genes file:', e)

    print('Wrote sampled least file:', sampled_path)
    print('Wrote audit:', audit_path)
    print('Sampled size:', len(sampled), 'target:', target_total)

    # Also (re)create the training data joined with features so that
    # `example/data_preprocessing/output/training_data_all_features.csv` and the
    # EDA input file are updated immediately after sampling. This mirrors the
    # join performed in the merge step so downstream ML steps see the filtered
    # least-likely genes without needing to re-run the merge script.
    try:
        training_df = pd.read_csv(config.training_genes, sep='\t')
        ml_df = pd.read_csv(config.all_genes_all_features_unprocessed)
        result_df = pd.merge(training_df, ml_df, on='Gene', how='left')
        # drop duplicates by Gene so training table has one row per gene
        if 'Gene' in result_df.columns:
            result_df.drop_duplicates(subset=['Gene'], inplace=True)
        else:
            result_df.drop_duplicates(inplace=True)
        # write both the ML features file and the EDA input file
        try:
            result_df.to_csv(config.training_genes_features, index=False)
            result_df.to_csv(config.training_data_all_features_eda, index=False)
            print('Wrote training features to:', config.training_genes_features)
            print('Wrote EDA training input to:', config.training_data_all_features_eda)
        except Exception as e:
            print('Warning: could not write joined training data files:', e)
    except Exception as e:
        print('Warning: could not join training_genes with all features:', e)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Stratified undersampling of least-likely genes')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--n_bins', type=int, default=3)
    parser.add_argument('--force', action='store_true', help='force sampling even if least <= combined')
    args = parser.parse_args()
    main(seed=args.seed, n_bins=args.n_bins, force=args.force)