process VALIDATE_INPUTS {
    tag "$params.dataset"
    label 'process_single'

    input:
    path protein_groups
    path design
    val  contrasts

    output:
    path "design.validated.csv", emit: design
    path "contrasts.txt",        emit: contrasts
    path "versions.yml",         emit: versions

    script:
    """
    export HOME="\$PWD"
    validate_inputs.R \\
        --protein_groups ${protein_groups} \\
        --design ${design} \\
        --contrasts "${contrasts}"
    """
}
