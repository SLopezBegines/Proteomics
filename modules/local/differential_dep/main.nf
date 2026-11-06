process DIFFERENTIAL_DEP {
    tag "$params.dataset"
    label 'process_medium'

    publishDir "${params.outdir}/${params.dataset}", mode: params.publish_dir_mode,
        saveAs: { fn -> fn.startsWith('out/') ? fn.substring(4) : null }

    input:
    path prev_state, stageAs: 'prev_state'

    output:
    path "state",                 emit: state
    path "out/**", type: 'file',  emit: published
    path "versions.yml",          emit: versions

    script:
    """
    export HOME="\$PWD"
    mkdir state && cp -rL ${prev_state}/. state/
    touch .marker && sleep 1
    differential_dep.R --state state --adjusted_method ${params.adjusted_method}

    mkdir out
    (cd state && find . -type f -newer ../.marker ! -name pipeline_context.rds -exec cp --parents -t ../out {} +)
    """
}
