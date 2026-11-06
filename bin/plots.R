#!/usr/bin/env Rscript
# Volcano / heatmap / PCA figures, Venn diagrams and summary statistics.
# Wraps code/05_Plots.R and code/17_summary_stats_proteomics.R (Venn diagrams: see note below).
this_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
source(file.path(this_dir, "pipeline_common.R"))

a <- parse_args(list(barplots = "0"))
if (is.null(a$state)) die("Missing --state")
a$state <- normalizePath(a$state, mustWork = TRUE)
ctx <- load_context(a$state)
setup_pipeline(ctx, a$state)
Exp_design <- ctx$exp_design
comparisons <- ctx$comparisons

print_barplots <- as.integer(a$barplots)
source("../code/05_Plots.R")

# NOTE: code/16_venn_diagram.R is intentionally NOT run. It filters on "<contrast>_diffexpressed",
# but code/04 now writes "<contrast>_diffexpressed_qval|pval|padj_holm_pval", so it fails on current
# results (the Rmd chunk fails the same way). Re-enable once that script is updated.

source("../code/17_summary_stats_proteomics.R")
data_results <- readxl::read_xlsx(file.path(output_path, "tables", "data_results.xlsx"))
list_results <- list(data_results)
names(list_results) <- ctx$dataset
all_stats <- apply_stats_to_list(list_results, comparisons)
print(all_stats)
writexl::write_xlsx(all_stats, file.path(output_path, "tables", "all_proteomic_stats.xlsx"))
# NOTE: plot_proteomic_stats() is not called: code/17 still uses save_plot(output_dir = ...), the
# pre-refactor signature (current save_plot() takes subdir). The statistics table above is still written.

save_context(ctx, a$state)
write_versions("PLOTS")
