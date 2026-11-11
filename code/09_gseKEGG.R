# gseKEGG: GSEA against KEGG pathways ####

source("../code/00_packages.R")

#' GSEA on KEGG pathways for the UP / DOWN sets of each comparison
#'
#' Ranks each set by fold change and runs clusterProfiler::gseKEGG(). Results are
#' saved to RData/ and tables/, lollipop plots to figures/gseKEGG.
#'
#' @param data_results       Table from run_de_analysis(); read from `output_path` if NULL.
#' @param comparisons        Character vector of comparison names.
#' @param diffexpressed_type "qval", "pval" or "padj_holm_pval".
#' @param top_n_FC           NULL = all DE proteins; n = the n with the largest |FC| per set.
#' @param kegg_organism      KEGG species code (params$kegg_organism, e.g. "mmu").
#' @param keyType            Gene ID type of the sets ("uniprot").
#' @param minGSSize,maxGSSize,pvalueCutoff,pAdjustMethod Passed to gseKEGG().
#' @param plot_lollipop      Draw and save the lollipop plots.
#' @param output_path        Project output folder.
#' @return List of gseaResult objects, one per set.
run_gseKEGG <- function(data_results = NULL,
                        comparisons,
                        diffexpressed_type = c("qval", "pval", "padj_holm_pval"),
                        top_n_FC = NULL,
                        kegg_organism = params$kegg_organism,
                        keyType = "uniprot",
                        minGSSize = 3,
                        maxGSSize = 800,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "none",
                        plot_lollipop = TRUE,
                        output_path = params$output_path) {
  diffexpressed_type <- match.arg(diffexpressed_type)
  if (is.null(data_results)) {
    data_results <- readxl::read_xlsx(file.path(output_path, "tables", "data_results.xlsx"))
  }

  sets <- prepare_DE_gene_sets(
    data_results = data_results, comparisons = comparisons,
    direction = c("UP", "DOWN"), independent_UPDOWN = TRUE,
    id_cols = c("name", "ID"), flatten = FALSE,
    diffexpressed_type = diffexpressed_type, top_n = top_n_FC, keep_stats = TRUE
  )
  gene_lists <- prepare_gse_gene_lists(sets)

  results_list_gseKEGG <- lapply(gene_lists, function(gl) {
    clusterProfiler::gseKEGG(
      geneList = gl, organism = kegg_organism, keyType = keyType,
      minGSSize = minGSSize, maxGSSize = maxGSSize,
      pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod
    )
  })
  save(results_list_gseKEGG, file = file.path(output_path, "RData", "results_list_gseKEGG_ALL.RData"))

  results_df <- purrr::map(results_list_gseKEGG, ~ .x@result)
  writexl::write_xlsx(results_df, file.path(output_path, "tables", "results_df_gseKEGG_ALL.xlsx"))

  if (plot_lollipop) {
    for (nm in names(results_df)) {
      tag <- gsub("^names_", "", nm)
      p <- plot_gse_lolliplot(results_df[[nm]], paste("Plot for", gsub("_", " ", tag)), x_label = "KEGG pathway")
      if (!is.null(p)) save_plot(paste0("KEGG_lolliplot_", tag), p, subdir = "gseKEGG", width = 7, height = 7)
    }
  }

  invisible(results_list_gseKEGG)
}
