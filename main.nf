#!/usr/bin/env nextflow
/*
 * Proteomics DEP pipeline: label-free LFQ proteomics (MaxQuant) with DEP/limma.
 * https://github.com/SLopezBegines/Proteomics  (branch docker-nextflow)
 */

nextflow.enable.dsl = 2

include { VALIDATE_INPUTS   } from './modules/local/validate_inputs/main'
include { PREPROCESS_QC     } from './modules/local/preprocess_qc/main'
include { DIFFERENTIAL_DEP  } from './modules/local/differential_dep/main'
include { PLOTS             } from './modules/local/plots/main'
include { SOFTWARE_VERSIONS } from './modules/local/software_versions/main'

// Fail early, with every problem listed at once, before any container is started.
def validateParams() {
    def errors = []
    if (!params.protein_groups) {
        errors << "--protein_groups is required (MaxQuant proteinGroups .txt/.tsv/.csv, or a pre-processed .xlsx)"
    } else if (!file(params.protein_groups).exists()) {
        def hint = params.protein_groups.toString().contains('testdata')
            ? " (the test subset is not versioned: see 'Test data' in the README)" : ""
        errors << "--protein_groups: file not found: ${params.protein_groups}${hint}"
    }
    if (!params.design) {
        errors << "--design is required (CSV with columns: label, condition, replicate[, batch])"
    } else if (!file(params.design).exists()) {
        errors << "--design: file not found: ${params.design}"
    }
    if (!(params.organism in ['human', 'mouse', 'zebrafish'])) {
        errors << "--organism must be one of human, mouse, zebrafish (got '${params.organism}')"
    }
    if (!(params.mnar_var in ['zero', 'MinProb', 'QRILC'])) {
        errors << "--mnar_var must be one of zero, MinProb, QRILC (got '${params.mnar_var}')"
    }
    if (!(params.adjusted_method in ['BH', 'BY', 'bonferroni', 'holm', 'fdr', 'none'])) {
        errors << "--adjusted_method must be one of BH, BY, bonferroni, holm, fdr, none (got '${params.adjusted_method}')"
    }
    if (!(params.fraction_na instanceof Number) || params.fraction_na < 0 || params.fraction_na > 1) {
        errors << "--fraction_na must be a number between 0 and 1 (got '${params.fraction_na}')"
    }
    if (!(params.p_val instanceof Number) || params.p_val <= 0 || params.p_val > 1) {
        errors << "--p_val must be a number in (0, 1] (got '${params.p_val}')"
    }
    if (!(params.p_val_low instanceof Number) || params.p_val_low <= 0 || params.p_val_low > 1) {
        errors << "--p_val_low must be a number in (0, 1] (got '${params.p_val_low}')"
    }
    if (!(params.fc instanceof Number) || params.fc < 0) {
        errors << "--fc (minimum |log2 fold change|) must be a number >= 0 (got '${params.fc}')"
    }
    if (params.contrasts && !(params.contrasts.toString().split(',').every { it.trim() ==~ /^.+_vs_.+$/ })) {
        errors << "--contrasts must be a comma-separated list of <condition1>_vs_<condition2> (got '${params.contrasts}')"
    }
    if (errors) {
        error "Invalid parameters:\n  - " + errors.join("\n  - ")
    }
}

workflow PROTEOMICS_DEP {
    main:
    ch_protein_groups = channel.fromPath(params.protein_groups, checkIfExists: true)
    ch_design         = channel.fromPath(params.design, checkIfExists: true)

    VALIDATE_INPUTS(ch_protein_groups, ch_design, params.contrasts ?: '')
    PREPROCESS_QC(ch_protein_groups, VALIDATE_INPUTS.out.design, VALIDATE_INPUTS.out.contrasts)
    DIFFERENTIAL_DEP(PREPROCESS_QC.out.state)
    PLOTS(DIFFERENTIAL_DEP.out.state)

    ch_versions = VALIDATE_INPUTS.out.versions
        .mix(PREPROCESS_QC.out.versions, DIFFERENTIAL_DEP.out.versions, PLOTS.out.versions)
        .collectFile(name: 'collated_versions.yml', newLine: true, sort: 'deep')
    SOFTWARE_VERSIONS(ch_versions)
}

workflow {
    validateParams()
    PROTEOMICS_DEP()
}
