
```
conda activate GenePrioritiser_env 
pip install pyensembl
pyensembl install --release 105 --species homo_sapiens
pip install --force-reinstall scikit-learn==1.4.2
pip install --force-reinstall scipy==1.11.4
pip install --force-reinstall numpy==1.23.0

```

Inputs/Requirements:
- GWAS data in format of GWAS catalog format summary statistics
    - Header: "MarkerName Allele1 Allele2 Freq1 Effect StdErr P TotalSampleSize N_effective"
- Ensure GWAS file name ends in phenotype as file name suffix (e.g. "GWAS_Evangelou_DBP.txt.gz")
- Include a text file list of most likely/known disease-specific genes and a list of probable disease-specific genes
    - Example: /GenePrioritiser/example/data_preprocessing/input/most_likely_genes.tsv
- Check database folder for disease-specific features (add new folders for new features) - all datasets need a 'Gene' column with HGNC gene symbols.
- Variant level data needs name format {filename}_{phenotype}.csv
    - It needs to go in folder: /GenePrioritiser/example/data_preprocessing/output/variants

Optional input:
'/GenePrioritiser/example/data_preprocessing/input/Additional_genes_to_filter.txt' can be made to further fitler least likely genes - identify genes via any requirements (e.g. genes with LD in the GWAS) for additional refining of least likely genes

# To Run:
```
nextflow run main.nf -profile conda
```