library(tidyverse)
library(magrittr)
library(data.table)
library(enrichR)
library(org.Hs.eg.db)

genes_all <- fread(least_likely_genes)
genes_all <- genes_all %>% dplyr::rename(Gene = Gene_Symbol)
genes_all <- dplyr::select(genes_all, Gene)
genes_all <- genes_all %>% distinct(Gene, .keep_all = TRUE)

db_list <- listEnrichrDbs()

db_list_df <- tibble::tibble(lib_names=db_list$libraryName)

dput(db_list_df$lib_names)

dbs <- c( "Chromosome_Location", "Human_Phenotype_Ontology", 
         "OMIM_Expanded",  
        "GO_Biological_Process_2023", "GO_Cellular_Component_2023", 
         "GO_Molecular_Function_2023",  "GWAS_Catalog_2023", "WikiPathway_2023_Human", 
         "KEGG_2021_Human",  "MGI_Mammalian_Phenotype_Level_4_2021", 
        "ClinVar_2019", 
         "DisGeNET", "Elsevier_Pathway_Collection")


res_all <- parallel::mclapply(genes_all$Gene, function(x) enrichr(x, dbs), mc.cores=20) %>% set_names(genes_all$Gene)

res_bind <- map(res_all, function(x) map(1:length(x), function(i) {
  print(names(x[i]))
  #map(x, function(i) print(names(i)))
  
  }))
  

res_bind <- map(1:length(res_list), function(x) {
  
  map(1:length(res_list[[x]]), function(i) {
  
  map(1:length(res_list[[x]][[i]]), function(j) {
    
    tibble(paste0(res_list[[x]][[i]][[j]]$Term, collapse=', ')) %>% set_names(names(res_list[[x]][[i]][j]))
    }) %>% bind_cols() }) %>% set_names(names(res_list[[x]])) %>% bind_rows(.id='Gene')})
  
  
res_bind_df <- res_bind %>% bind_rows()

#Add full gene name

gene_names_df <- select(org.Hs.eg.db, unique(res_bind_df$Gene), 
                        c("GENENAME"), "ALIAS") %>% distinct(ALIAS, .keep_all=TRUE) %>% dplyr::rename(`Gene name`=GENENAME)


res_bind_df_final <- res_bind_df %>% left_join(gene_names_df, by=c('Gene'='ALIAS')) %>% dplyr::select(`Gene symbol`=Gene, `Gene name`, everything())

count_bp_related <- res_bind_df_final %>%
    filter(str_detect(`Gene name`, regex("blood pressure|hypertension|cardio|cardiac|heart", ignore_case = TRUE))) %>%
    nrow()

print("Number of least likely genes enriched in BP or CVD pathways")
print(count_bp_related)

#fwrite(res_bind_df_final, paste0(enrichment_path, '/All_genes_annotated_with_EnrichR_', Sys.Date(), '.csv'))
