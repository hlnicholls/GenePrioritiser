# Gene Prioritisation Pipeline

Disease-specific gene prioritisation pipeline using GWAS summary statistics and user-defined gene-level features to prioritise genes likely involved in disease.

## Installation:

```
git clone https://github.com/hlnicholls/GenePrioritiser
cd GenePrioritiser

conda env create -f GenePrioritiser_env.yml
conda activate GenePrioritiser_env
python -m pip install --force-reinstall scikit-learn==1.4.2 scipy==1.11.4
pip install --force-reinstall numpy==1.23.0
conda install -c conda-forge scikit-optimize
conda install -c bioconda bedtools htslib
conda install -c conda-forge parallel
pip install pybedtools intervaltree requests

# R requirements: R (>= 4.0 recommended)
Rscript install_R_packages.R
```

## Requirements:

### Download GenePrioritiser-db
**Several databases need to be downloaded to provide gene annotation and model features to run the pipeline. These can be be downloaded from: [huggingface.co/datasets/hlnicholls/GenePrioritiser-db](https://huggingface.co/datasets/hlnicholls/GenePrioritiser-db)**

This dataset includes the `/databases` which needs to be downloaded to the repository root.

The gene annotation files need to be added to `/utils`

The dataset also includes an example which can be downloaded and run with this repository's current config settings (replace `/example` with the download `/example` from huggingface to run with the full data)

To download the huggingface dataset, you can use either Git LFS or the Hugging Face CLI/Python API.

#### Download postgwas-db with Git LFS
```bash
git lfs install --skip-repo
git clone https://huggingface.co/datasets/hlnicholls/GenePrioritiser-db
```

#### Or in Python:
```python
from huggingface_hub import snapshot_download
snapshot_download(
    repo_id="hlnicholls/GenePrioritiser-db",
    repo_type="dataset",
    local_dir="./databases"
)
```

### Inputs
- Add your GWAS data to ```/input/full_gwas```
- Update the ```/config/config.py``` file for your updated variables (e.g. GWAS file path```)
- GWAS data must be GRCh37-based harmonized summary statistics from GWAS Catalog (or similar)
    - File header/column names: ```chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value rsid n``` (or compatible variants)
    - Supported alternative column names are automatically detected:
      - **Chromosome**: `chromosome`, `chrom`, `chr`, `chromosome_name`
      - **Position**: `base_pair_location`, `basepair`, `bp`, `position`, `pos`, `genpos`
      - **Effect size**: `beta`, `effect`, `BETA`, `Effect`
      - **P-value**: `p_value`, `pvalue`, `p_val`, `pval`, `p.value`, `p`, `P`
- Ensure GWAS file name ends in ```_{phenotype}``` as file name suffix (e.g. "SBP_GCST90310294.h.tsv.gz" where SBP is the phenotype)
- Include a text file list of labelled most likely/known disease-specific genes and a list of probable disease-specific genes
    - Example (most likely genes are those validated by an expert clinician as interacting with BP drugs, probable genes are annotated by those that interact with any BP drugs in Open Targets/CHEMBL): 
        - ```/GenePrioritiser/example/data_preprocessing/input/most_likely_genes.tsv```
        - ```/GenePrioritiser/example/data_preprocessing/input/probable_genes.tsv```
        
- Check database folder and include disease-specific features (add new folders for new features) - all datasets need a ```Gene``` column with HGNC gene symbols.
- Any added variant level data for processing as a feature needs name format {filename}_{phenotype}.csv
    - It needs to go in folder: /GenePrioritiser/example/data_preprocessing/output/variants
    - The median value will be taken per gene

- Optional input can be made to further filter least likely genes identified in the pipeline for model training. This needs to be set to `least_likely_extra_filter` in the `config.py` (for example: `example/data_preprocessing/input/BP_loci_Apr2020_LDr2-8_500kb.csv`) - this can identify genes via any requirements (e.g. genes with LD in the GWAS) for additional refining of least likely genes. This will reduce your list of least likely genes based on their PPIs with the input additional gene list.

## Expected GWAS file column name format

The GWAS summary statistics should be GRCh37-based harmonized files from GWAS Catalog (or similar source). The following columns are **required**. The script automatically detects compatible column names using case-insensitive matching.

**Required columns:**

| Column Name | Alternative Names | Description |
|---|---|---|
| chromosome | `chrom`, `chr`, `chromosome_name` | Chromosome number (1-22, X, Y) |
| base_pair_location | `basepair`, `bp`, `position`, `pos`, `genpos` | Genomic position in GRCh37 |
| p_value | `pvalue`, `p_val`, `pval`, `p.value`, `p`, `P` | P-value for the association |
| beta | `effect`, `Effect`, `BETA` | Effect size (beta coefficient or log odds ratio) |

**Optional but recommended columns:**

| Column Name | Description |
|---|---|
| effect_allele | Alternative allele at the locus |
| other_allele | Reference allele at the locus |
| effect_allele_frequency | Frequency of the effect allele in the study population |
| standard_error | Standard error of the effect size |
| rsid | dbSNP rsID for the variant |
| n | Sample size for this variant |

**Example GWAS file (tab-separated, gzip-compressed):**

```
chromosome	base_pair_location	effect_allele	other_allele	beta	standard_error	effect_allele_frequency	p_value	rsid	n
1	758351	G	A	0.0568	0.0584	0.1169	0.3309	rs12238997	730165
1	787290	C	T	0.0242	0.076	0.0975	0.7499	rs116030099	730165
1	794707	C	T	0.1523	0.0854	0.0514	0.07446	rs148120343	730165
1	796338	C	T	0.0387	0.0517	0.1256	0.4538	rs58276399	730165
```

Include your GRCh37 GWAS file(s) in `GenePrioritiser/input/full_gwas/`.

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

# Create output results folder directory structure (can be named as desired)
bash src/data_preprocessing/Step0_directory_creation.sh "./results"

# rename all /example/ files in config. to /results
# upload GWAS file to /results/data_preprocessing/input
# upload your most likely/probable genes files to /results/data_preprocessing/input
# ensure input file names align with your config.py settings

nextflow run main.nf -profile conda
```

## Outputs
 - Results are found in output subfolders. E.g., ```/GenePrioritiser/results/machine_learning/multiclass/output```

## Example

An example GWAS summary statistics file and training genes files are included in hugginface dataset (prioritising genes for blood pressure GWASs). Interactive results are also shown on huggingface: [huggingface.co/spaces/hlnicholls/BP-GWAS-Prioritise](https://huggingface.co/spaces/hlnicholls/BP-GWAS-Prioritise)