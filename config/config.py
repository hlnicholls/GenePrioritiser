
# 1. Two files need to be provided by the user:
gwas_path = [
    "/Users/hannahnicholls/GitHub/GenePrioritiser/example/data_preprocessing/input/HF_HRC_GWAS_UKBB_EUR_HF.txt.gz"
] # https://www.ebi.ac.uk/gwas/studies/GCST007715
most_likely_gene_path = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/data_preprocessing/input/most_likely_genes.tsv"
# Example format of most likely file:
# Gene	label
# ABCC9	most likely
# ACBD4	most likely
# ACE	most likely

# Download drugs for phenotype of interest from Open Targets:
ot_phenotype_drugs = '/Users/hannahnicholls/GitHub/GenePrioritiser/databases/opentargets/EFO_0000537-known-drugs.tsv'

# 2. Check file paths for all intermediate files:
gwas_processed_path = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/data_preprocessing/output/HF_HRC_GWAS_UKBB_EUR_HF.txt.gz"

input_directory = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/data_preprocessing/input"
variant_output_directory = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/data_preprocessing/output/variants"
variant_database_output_directory = "/Users/hannahnicholls/GitHub/GenePrioritiser/databases/variant_level"

database_path = "/Users/hannahnicholls/GitHub/GenePrioritiser/databases"

database_string_path = "/Users/hannahnicholls/GitHub/GenePrioritiser/databases/stringdb/"

least_likely_gene_path = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/data_preprocessing/input/least_likely_genes.tsv"

annotated_gwas = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/data_preprocessing/output/variants/Annotated_GWAS_HF.csv"
gene_types = "/Users/hannahnicholls/GitHub/GenePrioritiser/utils/hg19Rel92_AllgeneTypes_0kb.txt"

# Optional for advanced least likely gene filtering:
# loci_ld_file = '/../Genes_with_LD.txt'

# Probable genes defined by OT gene-drug interactions (not in most likely gene group)
probable_gene_path = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/data_preprocessing/input/probable_genes.tsv"
# Example format of probable file:
# Gene	label
# AKT2	probable
# AMH	probable
# ANGPT2	probable

training_genes = '/Users/hannahnicholls/GitHub/GenePrioritiser/example/data_preprocessing/input/training_genes.txt'

all_genes_all_features_unprocessed = '/Users/hannahnicholls/GitHub/GenePrioritiser/example/data_preprocessing/output/all_genes_merged_all_data.csv'

training_genes_features = '/Users/hannahnicholls/GitHub/GenePrioritiser/example/data_preprocessing/output/training_data_all_features.csv'
training_data_all_features_eda = '/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/eda/input/training_data_all_features.csv'

genes_to_prioritise = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/data_preprocessing/output/gwas_genes_to_prioritise.csv"


# EDA
eda_script_path = '/Users/hannahnicholls/GitHub/GenePrioritiser/src/machine_learning/eda/scripts'

feature_missingness_plot = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/eda/output/feature_missingness.png"

correlation_pairs_09 = '/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/eda/output/correlation_pairs_09.csv'
correlation_matrix_09 = '/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/eda/output/correlation_matrix_09.csv'

cleaned_training_data = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/eda/output/cleaned_training_data.csv"

cleaned_training_data_eda = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/eda/output/cleaned_training_data.csv"

cleaned_training_data_ml = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/multiclass/input/cleaned_training_data.csv"

ml_eval_metrics = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/multiclass/output/model_evaluation_metrics.csv"

ml_output_path = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/multiclass/output"

boruta_shap_script = "/Users/hannahnicholls/GitHub/GenePrioritiser/src/machine_learning/multiclass/scripts"
boruta_shap_plot = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/multiclass/output/borutashap_feature_importance_plot_no_class_weights.png"

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

boruta_shap_plot_class_weighted = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/multiclass/output/borutashap_feature_importance_plot_class_weighted.png"

balanced_accuracy_per_fold = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/multiclass/output/per_fold_balanced_accuracies.csv"
ml_eval_metrics_class_weighted = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/multiclass/output/model_evaluation_metrics_class_weighted.csv"
best_voting_model_pkl_file = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/multiclass/output/best_voting_model_fitted.pkl"
best_model_pkl_file = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/multiclass/output/best_model_fitted.pkl"

best_model_predictions = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/multiclass/output/best_model_predictions_all_gwas_genes_with_probabilities.csv"

all_imputed_features ="/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/multiclass/output/all_genes_imputed_features.csv"

best_model_predictions_all_genes = "/Users/hannahnicholls/GitHub/GenePrioritiser/example/machine_learning/multiclass/output/best_model_predictions_all_genes_with_probabilities.csv"