import os
import pandas as pd
import fnmatch
import sys
CURRENT_DIR = os.getcwd()
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
sys.path.append(CONFIG_DIR)
import config

# Directory containing input files
input_directory = config.variant_output_directory 
# Merged beta(s) file path
output_file = os.path.join(config.variant_database_output_directory, "merged_gene_median_variant_measures.csv")

# Initialize a list to store DataFrames
processed_dfs = []

# Define the filename pattern
filename_pattern = "variant_data_*.csv"

# Iterate over each file in the directory
for filename in os.listdir(input_directory):
    if fnmatch.fnmatch(filename, filename_pattern):
        # Construct the full file path
        file_path = os.path.join(input_directory, filename)

        # Skip the output file if it exists in the input directory
        if file_path == output_file:
            print(f"Skipping the existing merged file: {file_path}")
            continue

        # Read the CSV file
        df = pd.read_csv(file_path)

        # Trim whitespace from all column names
        df.columns = df.columns.str.strip()

        # Check for 'Gene' column and one other column (assume it's the measure column)
        gene_columns = [col for col in df.columns if col.lower() == 'gene']
        other_columns = [col for col in df.columns if col.lower() != 'gene']

        # Validate that there is exactly one 'Gene' column and one measure column
        if len(gene_columns) != 1 or len(other_columns) != 1:
            print(f"The data in {filename} needs to contain 1 'Gene' column and 1 variant measure column only to work.")
            continue  # Skip processing this file if the condition is not met

        # Identify the measure column
        gene_col = gene_columns[0]
        measure_col = other_columns[0].strip()

        # Clean whitespace and invalid characters from 'Gene' values
        df[gene_col] = df[gene_col].astype(str).str.strip().replace(r'^\s*$', pd.NA, regex=True)
        df = df.dropna(subset=[gene_col])  # Remove rows with NaN 'Gene'
        df = df[df[gene_col] != '']  # Remove rows with empty 'Gene'

        # Extract the first gene if multiple genes are listed
        df[gene_col] = df[gene_col].str.split().str[0]

        # Extract phenotype name by taking the substring after the last underscore and before '.csv'
        phenotype = filename.rsplit('_', 1)[-1].replace('.csv', '')

        # Rename the measure column to include the phenotype as a suffix
        new_measure_col = f"{measure_col}_{phenotype}_Median"
        df.rename(columns={measure_col: new_measure_col}, inplace=True)

        # Calculate median for each unique Gene
        median_df = df.groupby(gene_col, as_index=False)[new_measure_col].median()

        # Add the processed DataFrame to the list
        processed_dfs.append(median_df)

# Merge all processed DataFrames on 'Gene'
if processed_dfs:
    merged_df = processed_dfs[0]
    for df in processed_dfs[1:]:
        merged_df = pd.merge(merged_df, df, on='Gene', how='outer')

    # Ensure there is no extraneous whitespace, commas, or empty entries in 'Gene' column
    merged_df['Gene'] = merged_df['Gene'].astype(str).str.strip().replace(r'^\s*$', pd.NA, regex=True)
    merged_df = merged_df.dropna(subset=['Gene'])  # Remove rows with NaN 'Gene'
    merged_df = merged_df[merged_df['Gene'] != '']  # Remove rows with empty 'Gene'
    merged_df['Gene'] = merged_df['Gene'].str.replace(r'[^\w\s]', '', regex=True)  # Remove special characters

    # Save the merged output
    merged_df.to_csv(output_file, index=False)
    print("All files processed and merged. Output saved to:", output_file)
else:
    print("No valid files were processed.")
