import os

PROJECT_DIR = os.environ.get("projectDir", ".")

# ── Resource limits ──────────────────────────────────────────────────────────
# Maximum CPU cores and RAM available to the ML training scripts.
# Adjust these to match your server's available resources.
MAX_CORES = 16
MAX_RAM_GB = 64

# 1. GWAS input files (GRCh37-based, from GWAS Catalog or similar)
# These are harmonized summary statistics with columns:
#   chromosome, base_pair_location, beta, p_value (column names are auto-detected)
gwas_path = [
    os.path.join(PROJECT_DIR, "input/full_gwas/SBP_GCST90310294.tsv.gz"),
    os.path.join(PROJECT_DIR, "input/full_gwas/DBP_GCST90310295.tsv.gz"),
    os.path.join(PROJECT_DIR, "input/full_gwas/PP_GCST90310296.tsv.gz")
]

# Training gene lists (phenotype-specific known/probable/least-likely genes)
most_likely_gene_path = os.path.join(PROJECT_DIR, "input/training/most_likely_genes.tsv")
probable_gene_path = os.path.join(PROJECT_DIR, "input/training/probable_genes.tsv")
least_likely_gene_path = os.path.join(PROJECT_DIR, "input/training/least_likely_genes.tsv")

# Download drugs for phenotype of interest from Open Targets:
ot_phenotype_drugs = os.path.join(PROJECT_DIR, "databases/opentargets/EFO_0000537-known-drugs.tsv")

# 2. Output and intermediate file paths
input_directory = os.path.join(PROJECT_DIR, "results/data_preprocessing/output")
variant_output_directory = os.path.join(PROJECT_DIR, "results/data_preprocessing/output/variants")
variant_database_output_directory = os.path.join(PROJECT_DIR, "databases/variant_level")

database_path = os.path.join(PROJECT_DIR, "databases")
database_string_path = os.path.join(PROJECT_DIR, "databases/stringdb/")

# Gene annotation reference
gene_types = os.path.join(PROJECT_DIR, "utils/hg19Rel92_AllgeneTypes_0kb.txt")

# Probable gene selection options
# p-value threshold for declaring a SNP "significant" when selecting probable genes.
# Can be overridden by setting `probable_gene_pvalue_threshold` in this config.
probable_gene_pvalue_threshold = 0.01 #0.00001
# Control whether a gene must have a significant SNP in EVERY annotated GWAS file
# (intersection) or in ANY file (union). Default is union (False) which means a
# gene only needs at least one significant SNP in any file to qualify.
probable_gene_intersection = True
# Least-likely gene selection p-value threshold (default: 0.01)
least_likely_pvalue_threshold = 0.05
# Least-likely PPI exclusion mode:
#   "direct"               -> exclude only direct interactors (1 hop)
#   "direct_and_secondary" -> exclude direct + second-degree interactors (2 hops)
least_likely_ppi_mode = "direct_and_secondary"

training_genes = os.path.join(PROJECT_DIR, "results/data_preprocessing/input/training_genes.txt")

all_genes_all_features_unprocessed = os.path.join(PROJECT_DIR, "results/data_preprocessing/output/all_genes_merged_all_data.csv")

training_genes_features = os.path.join(PROJECT_DIR, "results/data_preprocessing/output/training_data_all_features.csv")
training_data_all_features_eda = os.path.join(PROJECT_DIR, "results/machine_learning/eda/input/training_data_all_features.csv")

genes_to_prioritise = os.path.join(PROJECT_DIR, "results/data_preprocessing/output/gwas_genes_to_prioritise.csv")

# EDA
eda_script_path = os.path.join(PROJECT_DIR, "src/machine_learning/eda")
feature_missingness_plot = os.path.join(PROJECT_DIR, "results/machine_learning/eda/output/feature_missingness.png")
correlation_pairs_09 = os.path.join(PROJECT_DIR, "results/machine_learning/eda/output/correlation_pairs_09.csv")
correlation_matrix_09 = os.path.join(PROJECT_DIR, "results/machine_learning/eda/output/correlation_matrix_09.csv")
cleaned_training_data = os.path.join(PROJECT_DIR, "results/machine_learning/eda/output/cleaned_training_data.csv")
cleaned_training_data_eda = os.path.join(PROJECT_DIR, "results/machine_learning/eda/output/cleaned_training_data.csv")

# Machine learning
cleaned_training_data_ml = os.path.join(PROJECT_DIR, "results/machine_learning/multiclass/input/cleaned_training_data.csv")
ml_eval_metrics = os.path.join(PROJECT_DIR, "results/machine_learning/multiclass/output/model_evaluation_metrics.csv")
ml_output_path = os.path.join(PROJECT_DIR, "results/machine_learning/multiclass/output")
boruta_shap_script = os.path.join(PROJECT_DIR, "src/machine_learning/multiclass")
boruta_shap_plot = os.path.join(PROJECT_DIR, "results/machine_learning/multiclass/output/borutashap_feature_importance_plot_no_class_weights.png")

xgb_parameters = {'max_depth': (1, 4), 'learning_rate': (0.01, 0.2, 'log-uniform'), 'n_estimators': (10, 50)}
lgbm_parameters = {"max_depth": (1, 4), "learning_rate": (0.01, 0.2, "log-uniform"), "n_estimators": (10, 50)}
catboost_parameters = {"iterations": (10, 50), 'learning_rate': (0.01, 0.2, 'log-uniform'), 'depth': (1, 4)}
rf_parameters = {'n_estimators': (10, 50), 'max_depth': (1, 4)}
lr_parameters = {'C': (1e-4, 1e2, 'log-uniform')}
svc_parameters = {'C': (1e-4, 1e2, 'log-uniform')}
nn_parameters = {
    'epochs': (10, 100),  # Number of epochs
    'batch_size': (16, 64),  # Batch size
}

boruta_shap_plot_class_weighted = os.path.join(PROJECT_DIR, "results/machine_learning/multiclass/output/borutashap_feature_importance_plot_class_weighted.png")
balanced_accuracy_per_fold = os.path.join(PROJECT_DIR, "results/machine_learning/multiclass/output/per_fold_balanced_accuracies.csv")
ml_eval_metrics_class_weighted = os.path.join(PROJECT_DIR, "results/machine_learning/multiclass/output/model_evaluation_metrics_class_weighted.csv")
best_voting_model_pkl_file = os.path.join(PROJECT_DIR, "results/machine_learning/multiclass/output/best_voting_model_fitted.pkl")
best_model_pkl_file = os.path.join(PROJECT_DIR, "results/machine_learning/multiclass/output/best_model_fitted.pkl")
best_model_predictions = os.path.join(PROJECT_DIR, "results/machine_learning/multiclass/output/best_model_predictions_all_gwas_genes_with_probabilities.csv")
all_imputed_features = os.path.join(PROJECT_DIR, "results/machine_learning/multiclass/output/all_genes_imputed_features.csv")
best_model_predictions_all_genes = os.path.join(PROJECT_DIR, "results/machine_learning/multiclass/output/best_model_predictions_all_genes_with_probabilities.csv")
