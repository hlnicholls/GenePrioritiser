nextflow.enable.conda = true

process step0_combine_multiple_GWAS {
    conda 'GenePrioritiser_env'
    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step0_combine_multiple_GWAS.py
    """
}

process step1_annotate_genes {
    conda 'GenePrioritiser_env'
    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step1_annotate_genes.py
    """
}

process step2_process_variant_level_data {
    conda 'GenePrioritiser_env'
    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step2_process_variant_level_data.py
    """
}

process step3_least_likely_gene_selection {
    conda 'GenePrioritiser_env'
    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step3_least_likely_gene_selection.py
    """
}

process step6_identify_training_genes {
    conda 'GenePrioritiser_env'
    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step6_identify_training_genes.py
    """
}

process step7_merge_all_databases_and_get_training_data {
    conda 'GenePrioritiser_env'
    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step7_merge_all_databases_and_get_training_data.py
    """
}

process step8_subset_genes_to_prioritise {
    conda 'GenePrioritiser_env'
    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step8_subset_genes_to_prioritise.py
    """
}

process eda_training_data {
    conda 'GenePrioritiser_env'
    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/machine_learning/eda/scripts/eda_training_data.py
    """
}

process model_benchmark {
    conda 'GenePrioritiser_env'
    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/machine_learning/multiclass/scripts/model_benchmark.py
    """
}

process model_class_weights_benchmark {
    conda 'GenePrioritiser_env'
    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/machine_learning/multiclass/scripts/model_class_weights_benchmark.py
    """
}

process best_model_prioritisation {
    conda 'GenePrioritiser_env'
    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/machine_learning/multiclass/scripts/best_model_prioritisation.py
    """
}

workflow {
    step0_combine_multiple_GWAS()
        .then(step1_annotate_genes)
        .then(step2_process_variant_level_data)
        .then(step3_least_likely_gene_selection)
        .then(step6_identify_training_genes)
        .then(step7_merge_all_databases_and_get_training_data)
        .then(step8_subset_genes_to_prioritise)
        .then(eda_training_data)
        .then(model_benchmark)
        .then(model_class_weights_benchmark)
        .then(best_model_prioritisation)
}
