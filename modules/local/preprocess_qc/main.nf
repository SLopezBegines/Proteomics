process PREPROCESS_QC {
    tag "$params.dataset"
    label 'process_high_memory'

    publishDir "${params.outdir}/${params.dataset}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.startsWith('out/') ? fn.substring(4) : null }

    input:
    path protein_groups
    path design
    path contrasts

    output:
    path "state",                 emit: state
    path "out/**", type: 'file',  emit: published
    path "versions.yml",          emit: versions

    script:
    """
    export HOME="\$PWD"
    mkdir state
    preprocess_qc.R \\
        --protein_groups ${protein_groups} \\
        --design ${design} \\
        --contrasts ${contrasts} \\
        --dataset ${params.dataset} \\
        --organism ${params.organism} \\
        --state state \\
        --fraction_na ${params.fraction_na} \\
        --factor_sd_impute ${params.factor_sd_impute} \\
        --factor_sd_scale ${params.factor_sd_scale} \\
        --mnar_var ${params.mnar_var} \\
        --p_val ${params.p_val} \\
        --p_val_low ${params.p_val_low} \\
        --fc ${params.fc}

    mkdir out
    (cd state && find . -type f ! -name pipeline_context.rds -exec cp --parents -t ../out {} +)
    """
}
