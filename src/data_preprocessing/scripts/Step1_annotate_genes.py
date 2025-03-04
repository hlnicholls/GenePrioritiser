import os
import gzip
import pandas as pd
from pyensembl import EnsemblRelease  # PyEnsembl is used to annotate genes. You may need to install this package.
import sys
CURRENT_DIR = os.getcwd()
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR))
CONFIG_DIR = os.path.join(ROOT_DIR, "config")
sys.path.append(CONFIG_DIR)
import config


# Initialize the EnsemblRelease (modify release number as needed for reference)
data = EnsemblRelease(105) # GRCh37

# Annotate the nearest gene within 10kb
def annotate_nearest_gene(chrom, pos):
    genes = data.genes_at_locus(contig=chrom, position=pos)
    # Filter genes by distance within 10kb
    nearby_genes = [gene for gene in genes if abs(gene.start - pos) <= 5000 or abs(gene.end - pos) <= 5000]
    if nearby_genes:
        # Return the first gene found within the range
        return nearby_genes[0].gene_name
    else:
        return 'NoGene'

# Process each file in the input directory
for filename in config.gwas_path: #os.listdir(config.input_directory):
    if filename.endswith(".txt.gz"):
        # Construct input file path
        input_file = os.path.join(config.input_directory, filename)
        
        # Extract phenotype name from filename
        phenotype = filename.split("_")[-1].replace(".txt.gz", "")
        
        # Construct output file paths
        output_file = os.path.join(config.variant_output_directory, f"Annotated_GWAS_{phenotype}.csv")
        variant_output_file = os.path.join(config.variant_output_directory, f"variant_data_{phenotype}.csv")
        
        # Read the gzipped file into a DataFrame
        with gzip.open(input_file, 'rt') as f:
            df = pd.read_csv(f, sep='\s+')
        
        print('Phenotype')
        print(df.columns)
        print('Input data:')
        print(df.head(5))
        df = df.rename(columns={
            "effect_allele": "Allele1",
            "EA": "Allele1",
            "A1": "Allele1",
            "noneffect_allele": "Allele2",
            "NEA": "Allele2",
            "A2": "Allele2",
            "chr": "CHROM",
            "Chromosome": "CHROM",
            "pos": "GENPOS",
            "position": "GENPOS",
            "BETA": "Effect",
            "beta": "Effect",
            "Beta": "Effect",
            "Pvalue" : "P",
            "pvalue" : "P",
            "p.value" : "P",
        })
        
        # Ensure 'MarkerName' is present, with case insensitive check
        marker_name_col = [col for col in df.columns if 'MarkerName' in col]
        if not marker_name_col:
            print(f"No 'MarkerName' column found in {filename}. Skipping.")
            continue
        marker_name_col = marker_name_col[0]  # Use the first matching column

        # Create the SNP column
        if not df[marker_name_col].str.contains("rs").any():
            df['SNP'] = df[marker_name_col].str.replace(':SNP', '', regex=False) + ':' + df['Allele1'].str.upper() + ':' + df['Allele2'].str.upper()
             # Extract CHROM and GENPOS from MarkerName
            df['CHROM'] = df[marker_name_col].str.split(':').str[0]
            df['GENPOS'] = df[marker_name_col].str.split(':').str[1].astype(int)
        else:
            df['SNP'] = df['CHROM'].astype(str) + ':' + df['GENPOS'].astype(str) + ':' + df['Allele1'].astype(str) + ':' + df['Allele2'].astype(str)

        
        print(f'Processing {filename}...')

        # Annotate genes with nearest gene within 10kb
        df['Gene'] = df.apply(lambda row: annotate_nearest_gene(row['CHROM'], row['GENPOS']), axis=1)
        df = df[df['Gene'] != 'NoGene']
        
        # Save the annotated data to a CSV file
        df.to_csv(output_file, index=False)
        print(f"Annotation complete for {filename}. Output saved to: {output_file}")

        # Create a second DataFrame with only 'Gene' and 'Effect' columns and save it
        if 'Effect' in df.columns:
            variant_df = df[['Gene', 'Effect']]
            variant_df.to_csv(variant_output_file, index=False)
            print(f"Variant data saved to: {variant_output_file}")
        else:
            print(f"No 'Effect' column found in the input data for {filename}.")
