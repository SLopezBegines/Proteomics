process SOFTWARE_VERSIONS {
    label 'process_single'

    publishDir "${params.outdir}/${params.dataset}/pipeline_info", mode: params.publish_dir_mode

    input:
    path versions

    output:
    path "software_versions.yml", emit: yml

    script:
    """
    echo '"Workflow":' > software_versions.yml
    echo '    nextflow: ${workflow.nextflow.version}' >> software_versions.yml
    echo '    ${workflow.manifest.name}: ${workflow.manifest.version}' >> software_versions.yml
    echo '    container: ${params.container}' >> software_versions.yml
    cat ${versions} >> software_versions.yml
    """
}
