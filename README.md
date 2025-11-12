# Gene Prioritisation Pipeline

Disease-specific gene prioritisation pipeline using GWAS summary statistics and various gene-level features to prioritise genes likely involved in disease.


```
conda env create -f GenePrioritiser_env.yml
conda activate GenePrioritiser_env
pip install --force-reinstall scikit-learn==1.4.2
pip install --force-reinstall scipy==1.11.4
pip install --force-reinstall numpy==1.23.0
pip install scikit-optimize
```

Inputs/Requirements:
- Add your GWAS data to ```/src/data_preprocessing/input```
- Update the ```/config/config.py``` file for your updated variables (e.g. GWAS file name```)
- GWAS data in format of GWAS catalog format summary statistics in GRCh37
    - File header/column names: ```MarkerName Allele1 Allele2 Freq1 Effect StdErr P TotalSampleSize N_effective```
- Ensure GWAS file name ends in ```_{phenotype}``` as file name suffix (e.g. "GWAS_Evangelou_DBP.txt.gz" where DBP is the phenotype)
- Include a text file list of labelled most likely/known disease-specific genes and a list of probable disease-specific genes
    - Example (most likely genes are those validated by an expert clinician as interacting with BP drugs, probable genes are annotated by those that interact with any BP drugs in Open Targets/CHEMBL): 
        - ```/GenePrioritiser/example/data_preprocessing/input/most_likely_genes.tsv```
        - ```/GenePrioritiser/example/data_preprocessing/input/probable_genes.tsv```
        
- Check database folder and include disease-specific features (add new folders for new features) - all datasets need a ```Gene``` column with HGNC gene symbols.
- Any added variant level data for processing as a feature needs name format {filename}_{phenotype}.csv
    - It needs to go in folder: /GenePrioritiser/example/data_preprocessing/output/variants
    - The median value will be taken per gene

## Expected GWAS file column name format

The GWAS summary statistics should include the following columns (tab/whitespace-separated is fine). Example rows are shown below:

| MarkerName | Allele1 | Allele2 | Freq1 | Effect | StdErr | P | TotalSampleSize | N_effective |
|---|:---:|:---:|:---:|:---:|:---:|:---:|---:|---:|
| 10:100000625:SNP | a | g | 0.5663 | -0.0432 | 0.0174 | 0.01294 | 757601 | 756275 |
| 10:100000645:SNP | a | c | 0.7935 | 0.0263 | 0.0214 | 0.218 | 757599 | 754442 |
| 10:100001867:SNP | t | c | 0.0139 | -0.0424 | 0.0839 | 0.6128 | 742538 | 595603 |
| 10:100003242:SNP | t | g | 0.8835 | 0.0696 | 0.0269 | 0.009814 | 757600 | 752479 |

Include your GRCh37 GWAS file in `GenePrioritiser/src/data_preprocessing/input`.

## Expected training genes file format

The training genes file should be a two-column table (tab- or comma-separated) with the gene symbol and its label. Example:

| Gene | label |
|---|---|
| ABCC9 | most likely |
| ACBD4 | most likely |
| ACE | most likely |

Repeat for "probable" and "least likely" genes if curated separately.
- "least likely" genes can be annotated for you (determined by lack of PPIs with most likely/probable genes, and if provided SNPs not in LD with any target disease loci) during data preprocessing step.
- Include your training genes in ```/GenePrioritiser/src/data_preprocessing/input```

## Expected features file format in `databases` folder

Feature files placed under `GenePrioritiser/databases` should include a `Gene` column and the feature values. Example for a HIPred file:

| Gene | HIPred |
|---|---:|
| AGRN | 0.636561334 |
| NOC2L | 0.519413292 |
| B3GALT6 | 0.42872709 |
| C1orf159 | 0.146139458 |

Add new folders in `GenePrioritiser/databases` for any additional features.

Optional input:
'/GenePrioritiser/example/data_preprocessing/input/Additional_genes_to_filter.txt' can be made to further filter least likely genes for model training - identify genes via any requirements (e.g. genes with LD in the GWAS) for additional refining of least likely genes. This will reduce your list of least likely genes based on their PPIs with the input additional gene list.

# To Run:
```
conda activate GenePrioritiser_env 

mkdir -p results/data_preprocessing/{input,output} \
         results/machine_learning/eda/{input,output} \
         results/machine_learning/{multiclassinput,multiclass/output}

# rename all /example/ files in config. to /results
# upload GWAS file to /results/data_preprocessing/input
# upload your most likely/probable genes files to /results/data_preprocessing/input
# ensure input file names align with your config.py settings

nextflow run main.nf -profile conda
```

## Outputs
 - Results are found in output subfolders. E.g., ```/GenePrioritiser/results/machine_learning/multiclass/output```
