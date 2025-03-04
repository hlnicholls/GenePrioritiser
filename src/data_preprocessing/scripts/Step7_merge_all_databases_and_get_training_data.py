import os
import pandas as pd
import numpy as np
import sys
import os
CURRENT_DIR = os.getcwd()
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
sys.path.append(CONFIG_DIR)
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
        print(f"Error reading {file_path}: {e}")
        return None

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
        all_data_path = config.all_genes_all_features_unprocessed
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
        ml_df = pd.read_csv(ml_data_path)
        result_df = pd.merge(training_df, ml_df, on='Gene', how='left')
        result_df.drop_duplicates(inplace=True)
        for output_path in output_paths:
            result_df.to_csv(output_path, index=False)
        print(f"Joined data saved to {output_paths}")
    except Exception as e:
        print(f"Error during join operation: {e}")

run_data_merge()
