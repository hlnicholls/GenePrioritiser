process step1_annotate_genes {
    conda 'GenePrioritiser_env'
    // no input: Step 1 is the first preprocessing action (bedtools script)

    output:
    path "step1_output.txt"

    script:
    """
    export projectDir="${projectDir}"
    # Step 1 is a bash-driven bedtools annotation script
    bash "${projectDir}/src/data_preprocessing/Step1_annotate_genes_bedtools.sh" || true
    echo "step1_done" > step1_output.txt
    """
}

process step2_harmonise_genes {
    conda 'GenePrioritiser_env'
    input:
    path step1_out

    output:
    path "step2_output.txt"

    script:
    """
    export projectDir="${projectDir}"
    python "${projectDir}/src/data_preprocessing/Step2_harmonise_genes.py" > step2_output.txt
    """
}

process step3_process_variant_level_data {
    conda 'GenePrioritiser_env'
    input:
    path step2_out

    output:
    path "step3_output.txt"

    script:
    """
    export projectDir="${projectDir}"
    python "${projectDir}/src/data_preprocessing/Step3_process_variant_level_data.py" > step3_output.txt
    """
}

process step4_least_likely_gene_selection {
    conda 'GenePrioritiser_env'
    input:
    path step3_out

    output:
    path "step4_output.txt"

    script:
    """
    export projectDir="${projectDir}"
    python "${projectDir}/src/data_preprocessing/Step4_least_likely_gene_selection.py" > step4_output.txt
    """
}

process step5_geneset_enrichment {
    conda 'GenePrioritiser_env'
    input:
    path step4_out

    output:
    path "step5_output.txt"

    script:
    """
    export projectDir="${projectDir}"
    Rscript "${projectDir}/src/data_preprocessing/Step5_least_likely_geneset_enrichment.R" > step5_output.txt
    """
}

process step6_identify_training_genes {
    conda 'GenePrioritiser_env'
    input:
    path step5_out

    output:
    path "step6_output.txt"

    script:
    """
    export projectDir="${projectDir}"
    python "${projectDir}/src/data_preprocessing/Step6_identify_training_genes.py" > step6_output.txt
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
    export projectDir="${projectDir}"
    python "${projectDir}/src/data_preprocessing/Step7_merge_all_databases_and_get_training_data.py" > step7_output.txt
    """
}

process step8_downsample_least_likely_genes {
    conda 'GenePrioritiser_env'
    input:
    path step7_out

    output:
    path "step8_output.txt"

    script:
    """
    export projectDir="${projectDir}"
    python "${projectDir}/src/data_preprocessing/Step8_downsample_least_likely_genes.py" > step8_output.txt
    """
}

process step9_subset_genes_to_prioritise {
    conda 'GenePrioritiser_env'
    input:
    path step8_out

    output:
    path "step9_output.txt"

    script:
    """
    export projectDir="${projectDir}"
    python "${projectDir}/src/data_preprocessing/Step9_subset_genes_to_prioritise.py" > step9_output.txt
    """
}

process eda_training_data {
    conda 'GenePrioritiser_env'
    input:
    path step9_out

    output:
    path "eda_output.txt"

    script:
    """
    export projectDir="${projectDir}"
    python "\${projectDir}/src/machine_learning/eda/scripts/eda_training_data.py" > eda_output.txt
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
    export projectDir="${projectDir}"
    python "\${projectDir}/src/machine_learning/multiclass/scripts/model_benchmark.py" > benchmark_output.txt
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
    export projectDir="${projectDir}"
    python "\${projectDir}/src/machine_learning/multiclass/scripts/model_class_weights_benchmark.py" > weights_output.txt
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
    export projectDir="${projectDir}"
    python "\${projectDir}/src/machine_learning/multiclass/scripts/best_model_prioritisation.py" > final_output.txt
    """
}

workflow {
    // Data-preprocessing: run Step1 -> Step9 in order
    step1_out = step1_annotate_genes()
    step2_out = step2_harmonise_genes(step1_out)
    step3_out = step3_process_variant_level_data(step2_out)
    step4_out = step4_least_likely_gene_selection(step3_out)
    step5_out = step5_geneset_enrichment(step4_out)
    step6_out = step6_identify_training_genes(step5_out)
    step7_out = step7_merge_all_databases_and_get_training_data(step6_out)
    step8_out = step8_downsample_least_likely_genes(step7_out)
    step9_out = step9_subset_genes_to_prioritise(step8_out)

    // EDA and machine-learning stages
    eda_out = eda_training_data(step9_out)
    benchmark_out = model_benchmark(eda_out)
    weights_out = model_class_weights_benchmark(benchmark_out)
    best_model_prioritisation(weights_out)
}
