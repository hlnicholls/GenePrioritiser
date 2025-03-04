import pandas as pd
import gzip
import os
import sys
CURRENT_DIR = os.getcwd()
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR, "..", "..", ".."))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
sys.path.append(CONFIG_DIR)
import config


# GWAS file paths
input_files = config.gwas_path

# Output combined file
output_file = config.gwas_processed_path

# Create a list to store dataframes
dfs = []

for file in input_files:
    print('Appending file:', file)
    # Extract phenotype from the filename
    phenotype = os.path.basename(file).split("_")[-1].replace(".txt.gz", "")
    
    # Read the gzipped file
    with gzip.open(file, 'rt') as f:
        df = pd.read_csv(f, sep="\t")
    
    # Add the Phenotype column
    df["Phenotype"] = phenotype
    
    # Append to the list of dataframes
    dfs.append(df)

# Concatenate all dataframes
combined_df = pd.concat(dfs, ignore_index=True)

# Save the combined dataframe as a gzipped file
with gzip.open(output_file, 'wt') as f:
    combined_df.to_csv(f, sep="\t", index=False)

print(f"Combined file saved to {output_file}")
