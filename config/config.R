least_likely_genes <- "./example/data_preprocessing/input/sampling/least_likely_intermediate.tsv"

# Terms to consider as relevant disease terms when filtering enrichr pathway results
# These are matched case-insensitively against pathway/term names.
disease_terms <- c("blood pressure", "hypertension",  "hypotension", "cardio", "cardiac", "heart")

# Enrichr database list to query (can be extended)
enrichr_dbs <- c("GWAS_Catalog_2025", "WikiPathway_2023_Human", "KEGG_2021_Human", "DisGeNET", "GO_Biological_Process_2025")

# Output directory for enrichment results
enrichment_output_dir <- "./example/data_preprocessing/output/enrichr"

ml_plot_folder <- "./example/machine_learning/multiclass/output/"