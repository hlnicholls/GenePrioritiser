library(tidyverse)
library(magrittr)
library(data.table)
library(enrichR)
library(org.Hs.eg.db)

# Load project config (expects config/config.R to define least_likely_genes, disease_terms, enrichr_dbs, enrichment_output_dir)
config_path <- file.path('config', 'config.R')
if (!file.exists(config_path)) stop('config/config.R not found — please create it and set least_likely_genes and disease_terms')
source(config_path)

# Read least likely genes file (expecting a column named 'Gene' or single-column file)
genes_df <- fread(least_likely_genes)
if (!'Gene' %in% names(genes_df)) {
  # try common alternatives
  if ('Gene_Symbol' %in% names(genes_df)) genes_df <- genes_df %>% dplyr::rename(Gene = Gene_Symbol)
  else genes_df <- genes_df %>% dplyr::rename(Gene = names(genes_df)[1])
}
genes <- unique(genes_df$Gene)

message('Loaded ', length(genes), ' least-likely genes from: ', least_likely_genes)

# Decide which databases to query
dbs <- if (exists('enrichr_dbs')) enrichr_dbs else c('GWAS_Catalog_2025','WikiPathway_2023_Human','KEGG_2021_Human','DisGeNET','GO_Biological_Process_2025')

message('Querying Enrichr DBs and running per-gene enrichment (may take several minutes)...')

# Run enrichr per-gene using parallel workers (cross-platform)
# Use mclapply on Unix-like systems and a PSOCK cluster on Windows.
nc <- min(20L, parallel::detectCores(logical = TRUE))

# wrapper to call enrichr safely and return NULL on error
enrich_safe <- function(g) {
  tryCatch({
    # ensure enrichR is loaded in worker
    if (!requireNamespace('enrichR', quietly = TRUE)) stop('enrichR not available in worker')
    enrichR::enrichr(g, dbs)
  }, error = function(e) {
    message('enrichr error for ', g, ': ', conditionMessage(e))
    NULL
  })
}

if (.Platform$OS.type == 'windows') {
  cl <- parallel::makeCluster(nc)
  # ensure workers have enrichR available
  parallel::clusterEvalQ(cl, { library(enrichR); NULL })
  # export dbs to workers
  parallel::clusterExport(cl, varlist = c('dbs'), envir = environment())
  res_all <- parallel::parLapply(cl, genes, function(x) {
    tryCatch({ enrichR::enrichr(x, dbs) }, error = function(e) { message('worker error: ', conditionMessage(e)); NULL })
  })
  parallel::stopCluster(cl)
  names(res_all) <- genes
} else {
  res_all <- parallel::mclapply(genes, function(x) {
    enrich_safe(x)
  }, mc.cores = nc)
  names(res_all) <- genes
}

results_list <- list()
for (gene in names(res_all)) {
  gene_results <- res_all[[gene]]
  for (db_name in dbs) {
    if (!is.null(gene_results[[db_name]])) {
      db_table <- gene_results[[db_name]]
      # Keep only the standardized columns to avoid type mismatches when binding
      standardized_table <- db_table %>% transmute(
        Term = if ("Term" %in% colnames(db_table)) as.character(Term) else NA_character_,
        Genes = if ("Genes" %in% colnames(db_table)) as.character(Genes) else if ("Overlap" %in% colnames(db_table)) as.character(Overlap) else NA_character_,
        Database = db_name,
        Gene = gene
      )
      results_list[[length(results_list) + 1]] <- standardized_table
    }
  }
}

if (length(results_list) == 0) {
  message('No enrichr results returned. Exiting.')
  quit(status = 0)
}

# Combine and collapse terms per gene x database
results_df <- bind_rows(results_list)

# Sanity check: ensure expected columns are present and have correct names
expected_cols <- c('Gene', 'Database', 'Term', 'Genes')
missing_cols <- setdiff(expected_cols, names(results_df))
if (length(missing_cols) > 0) {
  stop('Enrichr results missing expected columns: ', paste(missing_cols, collapse = ', '))
}

# Group and collapse terms per gene x database (use group_by_at to avoid NSE issues)
results_collapsed <- results_df %>% group_by_at(vars('Gene', 'Database')) %>% summarise(Terms = paste(unique(na.omit(Term)), collapse = '; '), .groups = 'drop')

# Pivot to wide format: one row per gene with columns for each database containing concatenated terms
filtered_table <- results_collapsed %>% pivot_wider(names_from = Database, values_from = Terms, values_fill = NA)

# Safely rename the gene column to a human-friendly name; handle possible name variants
if ('Gene' %in% names(filtered_table)) {
  filtered_table <- filtered_table %>% dplyr::rename(`Gene symbol` = Gene)
} else if ('gene' %in% names(filtered_table)) {
  filtered_table <- filtered_table %>% dplyr::rename(`Gene symbol` = gene)
} else if ('Gene symbol' %in% names(filtered_table)) {
  # already present, do nothing
} else {
  stop('Expected a gene column in the pivoted results (tried Gene/gene/Gene symbol); found: ', paste(names(filtered_table), collapse = ', '))
}

# Build disease pattern and identify rows where any DB column contains a disease term
pattern <- paste(disease_terms, collapse = '|')
db_cols <- setdiff(names(filtered_table), 'Gene symbol')
match_mask <- apply(filtered_table[ , db_cols], 1, function(row) {
  any(sapply(row, function(cell) {
    if (is.na(cell) || cell == '') return(FALSE)
    stringr::str_detect(cell, regex(pattern, ignore_case = TRUE))
  }))
})

matched_genes <- filtered_table$`Gene symbol`[which(match_mask)]

if (!dir.exists(enrichment_output_dir)) dir.create(enrichment_output_dir, recursive = TRUE)
matched_out <- file.path(enrichment_output_dir, paste0('least_likely_genes_enriched_', Sys.Date(), '.tsv'))


fwrite(data.table(Gene = matched_genes), matched_out, sep='\t')
fwrite(data.table(Gene = setdiff(genes, matched_genes)), least_likely_genes, sep='\t')

message('Wrote matched genes to: ', matched_out)
message('Wrote filtered least-likely genes to: ', least_likely_genes)
message('Found ', length(matched_genes), ' least-likely genes that appear in disease-related enrichr terms')

