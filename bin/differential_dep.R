#!/usr/bin/env Rscript
# Differential expression with DEP/limma on the imputed SummarizedExperiment.
# Wraps code/04_data_analysis.R.
this_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
source(file.path(this_dir, "pipeline_common.R"))

a <- parse_args(list(adjusted_method = "BH"))
if (is.null(a$state)) die("Missing --state")
a$state <- normalizePath(a$state, mustWork = TRUE)
ctx <- load_context(a$state)
setup_pipeline(ctx, a$state)
Exp_design <- ctx$exp_design
comparisons <- ctx$comparisons

f <- file.path(output_path, "RData", "mixed_splited_imputation_prot_data.RData")
if (!file.exists(f)) die("Imputed object not found: ", f)
load(f)                                        # -> mixed_splited_imputation

source("../code/04_data_analysis.R")
data_analysis(mixed_splited_imputation, Exp_design, comparisons, output_path, adjusted_method = a$adjusted_method)

save_context(ctx, a$state)
write_versions("DIFFERENTIAL_DEP")
