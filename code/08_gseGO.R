# gseGO: GSEA against Gene Ontology ####

source("../code/00_packages.R")

#' GSEA on GO terms for the UP / DOWN sets of each comparison
#'
#' Ranks each set by fold change and runs clusterProfiler::gseGO(). Results are
#' saved to RData/ and tables/, plots to figures/gseGO.
#'
#' @param data_results       Table from run_de_analysis(); read from `output_path` if NULL.
#' @param comparisons        Character vector of comparison names.
#' @param diffexpressed_type "qval", "pval" or "padj_holm_pval".
#' @param top_n_FC           NULL = all DE proteins; n = the n with the largest |FC| per set.
#' @param ont                "ALL", "BP", "CC" or "MF".
#' @param organism           OrgDb package name (params$organism).
#' @param keyType            Gene ID type of the sets ("UNIPROT").
#' @param minGSSize,maxGSSize,pvalueCutoff,pAdjustMethod,scoreType Passed to gseGO().
#' @param plots              Any of "dot", "cnet", "heatmap", "ridge", "pmc", "lolli".
#' @param showCategory       Terms shown in dot plots.
#' @param output_path        Project output folder.
#' @return List of gseaResult objects, one per set.
run_gseGO <- function(data_results = NULL,
                      comparisons,
                      diffexpressed_type = c("qval", "pval", "padj_holm_pval"),
                      top_n_FC = NULL,
                      ont = "ALL",
                      organism = params$organism,
                      keyType = "UNIPROT",
                      minGSSize = 3,
                      maxGSSize = 800,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "none",
                      scoreType = "pos",
                      plots = c("dot", "cnet", "heatmap", "ridge", "pmc", "lolli"),
                      showCategory = 10,
                      output_path = params$output_path) {
  diffexpressed_type <- match.arg(diffexpressed_type)
  plots <- match.arg(plots, c("dot", "cnet", "heatmap", "ridge", "pmc", "lolli"), several.ok = TRUE)
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

  results_list_gseGO <- lapply(gene_lists, function(gl) {
    clusterProfiler::gseGO(
      geneList = gl, ont = ont, keyType = keyType, OrgDb = organism,
      minGSSize = minGSSize, maxGSSize = maxGSSize, pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod, scoreType = scoreType, verbose = TRUE
    )
  })
  save(results_list_gseGO, file = file.path(output_path, "RData", paste0("results_list_gseGO_", ont, ".RData")))

  results_df <- purrr::map(results_list_gseGO, ~ .x@result)
  writexl::write_xlsx(results_df, file.path(output_path, "tables", paste0("results_list_gseGO_", ont, ".xlsx")))

  # One builder per plot type; each returns a ggplot or NULL (no enriched terms)
  builders <- list(
    dot = list(label = "Dot plot", file = "GSE_dot", fun = function(x) {
      enrichplot::dotplot(x, showCategory = showCategory, split = ".sign") +
        facet_grid(. ~ .sign) + scale_color_viridis()
    }),
    cnet = list(label = "CNET plot", file = "GSE_cnet", fun = function(x) {
      enrichplot::cnetplot(clusterProfiler::simplify(x), foldChange = x@geneList)
    }),
    heatmap = list(label = "HeatMap plot", file = "GSE_heatmap", fun = function(x) {
      enrichplot::heatplot(x, foldChange = x@geneList)
    }),
    ridge = list(label = "Ridge plot", file = "GSE_ridge", fun = function(x) {
      enrichplot::ridgeplot(x) + scale_fill_viridis()
    }),
    pmc = list(label = "PMC plot", file = "GSE_PMC", fun = function(x) {
      # Number of PubMed articles mentioning the top 3 term descriptions
      enrichplot::pmcplot(x@result$Description[1:3], 2010:2026, proportion = FALSE)
    })
  )

  for (type in intersect(plots, names(builders))) {
    b <- builders[[type]]
    for (nm in names(results_list_gseGO)) {
      x <- results_list_gseGO[[nm]]
      if (nrow(x@result) == 0) {
        message("[gseGO] ", b$label, ": no enriched terms for ", nm)
        next
      }
      tag <- gsub("^names_", "", nm)
      p <- b$fun(x) + ggplot2::ggtitle(paste0(b$label, " for ", gsub("_", " ", tag)))
      save_plot(paste0(b$file, "_", tag), p, subdir = "gseGO", width = 7, height = 7)
    }
  }

  if ("lolli" %in% plots) {
    for (nm in names(results_df)) {
      tag <- gsub("^names_", "", nm)
      p <- plot_gse_lolliplot(results_df[[nm]], paste("Plot for", gsub("_", " ", tag)), x_label = "GO term")
      if (!is.null(p)) save_plot(paste0("GO_lolliplot_", tag), p, subdir = "gseGO", width = 7, height = 7)
    }
  }

  invisible(results_list_gseGO)
}
