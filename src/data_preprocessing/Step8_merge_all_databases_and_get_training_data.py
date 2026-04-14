import os
import pandas as pd
import numpy as np
import csv
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
if CONFIG_DIR not in sys.path:
    sys.path.insert(0, CONFIG_DIR)
import config

def read_file(file_path):
    try:
        if file_path.endswith('.txt') or file_path.endswith('.tsv'):
            with open(file_path, 'r') as f:
                first_line = f.readline()
                sep = '\t' if '\t' in first_line else ','
        else:
            sep = ','
        df = pd.read_csv(file_path, sep=sep, engine='python')
        df.columns = df.columns.str.replace('[^a-zA-Z0-9_]', '', regex=True)
        return df
    except Exception as e:
        # Retry with permissive parser settings for malformed quote characters.
        try:
            df = pd.read_csv(
                file_path,
                sep=sep,
                engine='python',
                on_bad_lines='skip',
                quoting=csv.QUOTE_NONE
            )
            df.columns = df.columns.str.replace('[^a-zA-Z0-9_]', '', regex=True)
            print(f"Warning: permissive parsing used for {file_path}")
            return df
        except Exception as e2:
            print(f"Error reading {file_path}: {e2}")
            return None


def collapse_duplicates_by_gene(df):
    """Collapse duplicate Gene rows by averaging numeric columns only."""
    if df.columns.duplicated().any():
        dup_names = df.columns[df.columns.duplicated()].tolist()
        print(f"Warning: duplicate column names found ({len(dup_names)}). Keeping first occurrence for: {sorted(set(dup_names))}")
        df = df.loc[:, ~df.columns.duplicated()]

    if 'Gene' not in df.columns:
        return df

    dup_count = df['Gene'].duplicated().sum()
    if dup_count <= 0:
        return df

    print(f"Found {dup_count} duplicate Gene entries in merged data - aggregating numeric features by mean to collapse duplicates.")

    # Keep a copy of Gene and coerce all feature columns to numeric; non-numeric
    # values (e.g. 'Y') become NaN and do not break mean aggregation.
    coerced = df.copy()
    feature_cols = [c for c in coerced.columns if c != 'Gene']
    for col in feature_cols:
        coerced[col] = pd.to_numeric(coerced[col], errors='coerce')

    numeric_cols = [c for c in feature_cols if pd.api.types.is_numeric_dtype(coerced[c])]
    if not numeric_cols:
        return coerced[['Gene']].drop_duplicates().reset_index(drop=True)

    collapsed = coerced.groupby('Gene', as_index=False)[numeric_cols].mean()
    return collapsed

def merge_data(base_dir, excluded_folders, file_types):
    main_df = None
    for root, dirs, files in os.walk(base_dir):
        dirs[:] = [d for d in dirs if d not in excluded_folders]
        for file in files:
            if file.endswith(file_types):
                file_path = os.path.join(root, file)
                df = read_file(file_path)
                if df is not None:
                    try:
                        if 'Gene' not in df.columns:
                            print(f"'Gene' column not found in {file_path}. Skipping this file.")
                            continue
                        if main_df is None:
                            main_df = df
                        else:
                            # Merge with suffixes to handle duplicate column names
                            main_df = pd.merge(main_df, df, how='outer', on='Gene', suffixes=('_left', '_right'))
                    except KeyError as e:
                        print(f"KeyError during merging with file: {file_path}. Error: {e}")
                        continue
    return main_df

def drop_non_numeric_columns(df):
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    # Ensure 'Gene' column is retained
    cols_to_keep = ['Gene'] + [col for col in numeric_cols if col != 'Gene']
    df = df[cols_to_keep]
    return df

def rename_columns(df, column_renames):
    return df.rename(columns=column_renames)

def filter_and_save(df, output_path):
    df.drop_duplicates(inplace=True)
    df.replace("", np.nan, inplace=True)
    df.to_csv(output_path, index=False)

def run_data_merge():
    base_dir = config.database_path
    excluded_folders = {"opentargets", "stringdb", "dgidb", "avana", "impc", "ipa_BP_2020"}
    file_types = ('.txt', '.csv', '.tsv')
    
    column_renames = {
        'RVIS': 'RVIS_genic_intolerance',
        'HIPred': 'HIPred_haploinsufficiency',
        'SDI': 'SDI_essentiality',
        'GDI_Score': 'gene_damage_index_score',
        'pLI': 'pLI_exac',
        'ui_panglaoDB': 'panglao_ubiquitousness_index',
        'ExomiserScore': 'ExomiserScore_increasedBP'
    }
    
    main_df = merge_data(base_dir, excluded_folders, file_types)
    if main_df is not None:
        main_df = rename_columns(main_df, column_renames)
        main_df = drop_non_numeric_columns(main_df)
        main_df = collapse_duplicates_by_gene(main_df)
        all_data_path = config.all_genes_all_features_unprocessed
        all_data_dir = os.path.dirname(all_data_path)
        if all_data_dir:
            os.makedirs(all_data_dir, exist_ok=True)
        filter_and_save(main_df, all_data_path)
        print(f"All merged databases saved to {all_data_path}")

        training_data_path = config.training_genes
        output_paths = [config.training_genes_features,
                        config.training_data_all_features_eda]
        join_training_data_with_features(training_data_path, all_data_path, output_paths)
    else:
        print("No data was merged. Please check the files and directories.")

def join_training_data_with_features(training_data_path, ml_data_path, output_paths):
    try:
        training_df = pd.read_csv(training_data_path, sep='\t')  # Assuming training data is tab-separated
        # Print counts per unique label in the training genes file (case-insensitive 'label')
        # This helps confirm how many 'most likely', 'probable', 'least likely' genes we have.
        label_col = None
        for col in training_df.columns:
            if col.lower() == 'label':
                label_col = col
                break
        if label_col is not None:
            counts = training_df[label_col].value_counts(dropna=False)
            print("Training genes counts per label:")
            for lbl, cnt in counts.items():
                print(f"  {lbl}: {cnt}")
        else:
            print(f"No 'label' column found in training data. Columns present: {list(training_df.columns)}")
        ml_df = pd.read_csv(ml_data_path)
        # If ml_df contains duplicate Gene rows, a left-merge will expand the
        # training set multiplicatively. Collapse duplicates in ml_df first by
        # aggregating numeric features (mean). If duplicates exist, warn the user.
        if 'Gene' in ml_df.columns:
            dup_ml = ml_df['Gene'].duplicated().sum()
            if dup_ml > 0:
                print(f"Warning: ML features file has {dup_ml} duplicate Gene rows - aggregating numeric features by mean to avoid duplicating training rows.")
                ml_df = collapse_duplicates_by_gene(ml_df)

        result_df = pd.merge(training_df, ml_df, on='Gene', how='left')
        # If training_df itself contained duplicate rows, drop exact duplicates
        result_df.drop_duplicates(inplace=True)
        for output_path in output_paths:
            out_dir = os.path.dirname(output_path)
            if out_dir:
                os.makedirs(out_dir, exist_ok=True)
            result_df.to_csv(output_path, index=False)
        print(f"Joined data saved to {output_paths}")
    except Exception as e:
        print(f"Error during join operation: {e}")

run_data_merge()
