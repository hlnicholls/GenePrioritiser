process step0_combine_multiple_GWAS {
    conda 'GenePrioritiser_env'
    output:
    path "step0_output.txt"

    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step0_combine_multiple_GWAS.py > step0_output.txt
    """
}

process step1_annotate_genes {
    conda 'GenePrioritiser_env'
    input:
    path step0_out
    
    output:
    path "step1_output.txt"

    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step1_annotate_genes.py > step1_output.txt
    """
}

process step2_process_variant_level_data {
    conda 'GenePrioritiser_env'
    input:
    path step1_out
    
    output:
    path "step2_output.txt"

    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step2_process_variant_level_data.py > step2_output.txt
    """
}

process step3_least_likely_gene_selection {
    conda 'GenePrioritiser_env'
    input:
    path step2_out
    
    output:
    path "step3_output.txt"

    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step3_least_likely_gene_selection.py > step3_output.txt
    """
}

process step6_identify_training_genes {
    conda 'GenePrioritiser_env'
    input:
    path step3_out
    
    output:
    path "step6_output.txt"

    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step6_identify_training_genes.py > step6_output.txt
    """
}

process step7_merge_all_databases_and_get_training_data {
    conda 'GenePrioritiser_env'
    input:
    path step6_out
    
    output:
    path "step7_output.txt"

    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step7_merge_all_databases_and_get_training_data.py > step7_output.txt
    """
}

process step8_subset_genes_to_prioritise {
    conda 'GenePrioritiser_env'
    input:
    path step7_out
    
    output:
    path "step8_output.txt"

    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/data_preprocessing/scripts/Step8_subset_genes_to_prioritise.py > step8_output.txt
    """
}

process eda_training_data {
    conda 'GenePrioritiser_env'
    input:
    path step8_out
    
    output:
    path "eda_output.txt"

    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/machine_learning/eda/scripts/eda_training_data.py > eda_output.txt
    """
}

process model_benchmark {
    conda 'GenePrioritiser_env'
    input:
    path eda_out
    
    output:
    path "benchmark_output.txt"

    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/machine_learning/multiclass/scripts/model_benchmark.py > benchmark_output.txt
    """
}

process model_class_weights_benchmark {
    conda 'GenePrioritiser_env'
    input:
    path benchmark_out
    
    output:
    path "weights_output.txt"

    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/machine_learning/multiclass/scripts/model_class_weights_benchmark.py > weights_output.txt
    """
}

process best_model_prioritisation {
    conda 'GenePrioritiser_env'
    input:
    path weights_out
    
    output:
    path "final_output.txt"

    script:
    """
    python /Users/hannahnicholls/GitHub/GenePrioritiser/src/machine_learning/multiclass/scripts/best_model_prioritisation.py > final_output.txt
    """
}

workflow {
    step0_out = step0_combine_multiple_GWAS()
    step1_out = step1_annotate_genes(step0_out)
    step2_out = step2_process_variant_level_data(step1_out)
    step3_out = step3_least_likely_gene_selection(step2_out)
    step6_out = step6_identify_training_genes(step3_out)
    step7_out = step7_merge_all_databases_and_get_training_data(step6_out)
    step8_out = step8_subset_genes_to_prioritise(step7_out)
    eda_out = eda_training_data(step8_out)
    benchmark_out = model_benchmark(eda_out)
    weights_out = model_class_weights_benchmark(benchmark_out)
    best_model_prioritisation(weights_out)
}
