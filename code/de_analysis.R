# de_analysis.R
#
# Differential expression analysis on an imputed proteomics
# SummarizedExperiment, using limma directly (no DEP layer).
#
# Replaces DEP::analyze_dep() + DEP::test_diff() + DEP::add_rejections() +
# DEP::get_results(), as called from data_analysis() in 04_data_analysis.R.
# Fits a ~0 + condition cell-means model, tests each "A_vs_B" string in
# `comparisons` as a limma contrast, and reproduces the SAME wide
# `data_results` column schema data_analysis() produced
# (<comparison>_ratio/_p.val/_p.adj/_q.val/_p.adj_holm/_diffexpressed_*/
# _dif_label_*, significance_*, significant, lfq_*, Score) so 05_Plots.R and
# the GO/KEGG/STRING scripts that read those columns by name keep working
# from the exported tables without changes.
#
# NOTE: unlike data_analysis(), this does NOT produce RData/dep_analysis.RData
# (there is no DEP SummarizedExperiment anymore). 05_Plots.R calls
# DEP::plot_pca()/plot_cor()/plot_heatmap()/plot_volcano()/plot_single()
# directly on that object, so it will break until migrated in a follow-up
# step — same pattern as when qc_proteomics.R replaced the older DEP-based
# QC/normalization/imputation scripts.
#
# Statistical helper functions are pure. run_de_analysis() is the
# orchestrator and the only function with side effects (Excel export, plot
# saving, assignment to .GlobalEnv) — same contract data_analysis() had.
#
# Function index:
#   build_de_design()
#   build_de_contrasts()
#   fit_de_model()
#   extract_contrast_results()
#   compute_qvalues()
#   build_results_table()
#   add_significance_labels()
#   plot_pvalue_distributions()
#   export_de_tables()
#   run_de_analysis()


#' Cell-means design matrix for a SummarizedExperiment
#'
#' ~0 + condition gives every condition its own coefficient (no reference
#' level), so manual "A_vs_B" contrasts map directly onto comparisons.
#'
#' @param se      SummarizedExperiment with a `condition` colData column
#' @param formula One-sided formula for model.matrix(); default ~0 + condition
#' @return model.matrix, columns named after the bare condition levels
#'   (the "condition" prefix model.matrix() adds is stripped)
build_de_design <- function(se, formula = ~0 + condition) {
  design <- stats::model.matrix(formula, data = as.data.frame(SummarizedExperiment::colData(se)))
  colnames(design) <- sub("^condition", "", colnames(design))
  design
}


#' Parse "A_vs_B" comparison strings into a limma contrasts matrix
#'
#' @param design      model.matrix from build_de_design()
#' @param comparisons Character vector like "CTRL_vs_WT" (split on the
#'   literal substring "_vs_")
#' @return limma contrasts matrix (condition levels x comparisons), columns
#'   named after `comparisons`
build_de_contrasts <- function(design, comparisons) {
  groups <- strsplit(comparisons, "_vs_", fixed = TRUE)
  bad <- lengths(groups) != 2
  if (any(bad)) {
    stop("Comparison(s) not in 'A_vs_B' form: ", paste(comparisons[bad], collapse = ", "))
  }
  missing_levels <- setdiff(unlist(groups), colnames(design))
  if (length(missing_levels) > 0) {
    stop("Condition(s) not found in design: ", paste(missing_levels, collapse = ", "))
  }

  contrast_exprs <- vapply(groups, function(g) paste(g[1], "-", g[2]), character(1))
  cm <- limma::makeContrasts(contrasts = contrast_exprs, levels = design)
  colnames(cm) <- comparisons
  cm
}


#' Fit protein-wise linear models + empirical Bayes moderation
#'
#' Operates on the "intensity" assay directly — it's already log2, fully
#' imputed (no NAs), continuous data, so no voom/weighting step is needed.
#'
#' @param se           SummarizedExperiment, assay "intensity"
#' @param design       model.matrix from build_de_design()
#' @param contrasts_mx contrasts matrix from build_de_contrasts()
#' @return eBayes-moderated MArrayLM fit, one coefficient per comparison
fit_de_model <- function(se, design, contrasts_mx) {
  mat <- SummarizedExperiment::assay(se, "intensity")
  fit <- limma::lmFit(mat, design)
  fit <- limma::contrasts.fit(fit, contrasts_mx)
  limma::eBayes(fit)
}


#' Per-protein logFC / p-value / BH-adjusted p-value for one contrast
#'
#' adj.P.Val is limma's own BH correction (topTable's default) — the
#' project's hardcoded primary FDR criterion (see CLAUDE.md).
#'
#' @param fit  eBayes fit from fit_de_model()
#' @param comp Comparison name (a column of the contrasts matrix)
#' @return data.frame(name, ratio, p.val, p.adj), one row per protein, in
#'   the same protein order as the SummarizedExperiment the fit came from
extract_contrast_results <- function(fit, comp) {
  tt <- limma::topTable(fit, coef = comp, number = Inf, sort.by = "none")
  data.frame(
    name  = rownames(tt),
    ratio = tt$logFC,
    p.val = tt$P.Value,
    p.adj = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
}


#' Storey q-values for one comparison's p-values
#'
#' @param pvals Numeric vector of raw p-values (may contain NA)
#' @return Numeric vector of q-values, same length/order as `pvals`
compute_qvalues <- function(pvals) {
  q <- rep(NA_real_, length(pvals))
  valid <- !is.na(pvals)
  q[valid] <- qvalue::qvalue(pvals[valid])$qvalues
  q
}


#' Assemble the wide per-comparison results table
#'
#' One row per protein. For each comparison in `comparisons`, adds
#' <comparison>_ratio, _p.val, _p.adj (limma BH), _q.val (Storey), and
#' _p.adj_holm (generic p.adjust(), method = adjusted_method — the "_holm"
#' suffix is kept for backward compatibility with downstream scripts even
#' though the method is configurable, matching the legacy schema).
#'
#' @param se              SummarizedExperiment the model was fit on (for name/ID/genename/intensity)
#' @param fit              eBayes fit from fit_de_model()
#' @param comparisons      Character vector of comparison names
#' @param adjusted_method  Method passed to p.adjust() for the `_p.adj_holm` column
#' @return Wide data.frame, one row per protein
build_results_table <- function(se, fit, comparisons, adjusted_method) {
  base <- data.frame(
    name     = rownames(se),
    ID       = SummarizedExperiment::rowData(se)$ID,
    genename = SummarizedExperiment::rowData(se)$genename,
    stringsAsFactors = FALSE
  )

  per_comparison <- lapply(comparisons, function(comp) {
    res <- extract_contrast_results(fit, comp)
    res$q.val <- compute_qvalues(res$p.val)
    res$p.adj_holm <- stats::p.adjust(res$p.val, method = adjusted_method)
    res <- res[, c("ratio", "p.val", "p.adj", "q.val", "p.adj_holm")]
    colnames(res) <- paste(comp, colnames(res), sep = "_")
    res
  })

  results <- do.call(cbind, c(list(base), per_comparison))

  lfq <- as.data.frame(SummarizedExperiment::assay(se, "intensity"))
  colnames(lfq) <- paste0("lfq_", colnames(lfq))

  cbind(results, lfq)
}


#' Add per-comparison UP/DOWN/NO labels, volcano labels, and significance flags
#'
#' For each of three adjustment flavours (raw p, Storey q, adjusted_method),
#' a protein is UP if ratio > lfc & stat < alpha, DOWN if ratio < -lfc &
#' stat < alpha. `significance_pval`/`_qval`/`_holmval` is TRUE if ANY
#' comparison is UP or DOWN under that flavour.
#'
#' `significant` replaces DEP::add_rejections(): TRUE if ANY comparison
#' passes the BH-adjusted p-value + lfc threshold. No _diffexpressed_adj/
#' _dif_label_adj columns are added for it, since the legacy schema never
#' had them either — it was only ever exposed as the `significant` flag.
#'
#' @param results     Wide table from build_results_table()
#' @param comparisons Character vector of comparison names
#' @param alpha       Significance threshold (p / q / adjusted-p < alpha)
#' @param lfc         |log2FC| threshold
#' @return `results` with added _diffexpressed_*, _dif_label_*, _Score,
#'   significant, and significance_* columns
add_significance_labels <- function(results, comparisons, alpha, lfc) {
  add_flavour <- function(results, comp, stat_col, diff_suffix, label_suffix) {
    diff_col  <- paste0(comp, "_diffexpressed_", diff_suffix)
    label_col <- paste0(comp, "_dif_label_", label_suffix)
    ratio <- results[[paste0(comp, "_ratio")]]
    stat  <- results[[stat_col]]

    call <- rep("NO", nrow(results))
    call[ratio > lfc & stat < alpha]  <- "UP"
    call[ratio < -lfc & stat < alpha] <- "DOWN"

    results[[diff_col]]  <- call
    results[[label_col]] <- ifelse(call != "NO", results$name, NA_character_)
    results
  }

  for (comp in comparisons) {
    results <- add_flavour(results, comp, paste0(comp, "_p.val"), "pval", "pval")
    results <- add_flavour(results, comp, paste0(comp, "_q.val"), "qval", "qval")
    results <- add_flavour(results, comp, paste0(comp, "_p.adj_holm"), "padj_holm_pval", "holmval")

    results[[paste0(comp, "_Score")]] <-
      results[[paste0(comp, "_ratio")]] * -log10(results[[paste0(comp, "_p.val")]])
  }

  any_flavour <- function(diff_suffix) {
    cols <- results[paste0(comparisons, "_diffexpressed_", diff_suffix)]
    apply(cols == "UP" | cols == "DOWN", 1, any)
  }
  results$significance_pval    <- any_flavour("pval")
  results$significance_qval    <- any_flavour("qval")
  results$significance_holmval <- any_flavour("padj_holm_pval")

  sig_adj <- sapply(comparisons, function(comp) {
    ratio <- results[[paste0(comp, "_ratio")]]
    padj  <- results[[paste0(comp, "_p.adj")]]
    (ratio > lfc | ratio < -lfc) & padj < alpha
  })
  results$significant <- apply(sig_adj, 1, any, na.rm = TRUE)

  results
}


#' Faceted histograms of p / adjusted-p / q / adjusted_method distributions
#'
#' @param results         Wide table from add_significance_labels()
#' @param adjusted_method Label for the p.adj_holm panel title
#' @return patchwork object (4 stacked histogram panels)
plot_pvalue_distributions <- function(results, adjusted_method) {
  histogram_panel <- function(suffix_regex, suffix, xlab, title) {
    cols <- grep(suffix_regex, names(results), value = TRUE)
    df <- results[cols] %>%
      tidyr::pivot_longer(cols = dplyr::everything(), names_to = "comparison", values_to = "value") %>%
      dplyr::mutate(comparison = stringr::str_remove(comparison, suffix))

    ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(bins = 100, fill = "steelblue", color = "black") +
      ggplot2::facet_wrap(~comparison, scales = "free_y") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = title, x = xlab, y = "Frequency")
  }

  p1 <- histogram_panel("_p\\.val$", "_p.val$", "P-value", "P-value distribution")
  p2 <- histogram_panel("_p\\.adj$", "_p.adj$", "Adjusted p-value", "Adjusted p-value distribution (BH)")
  p3 <- histogram_panel("_q\\.val$", "_q.val$", "Q-value", "Q-value distribution (Storey)")
  p4 <- histogram_panel(
    "_p\\.adj_holm$", "_p.adj_holm$",
    paste0("Adjusted p-value (", adjusted_method, ")"),
    paste0("Adjusted p-value distribution (", adjusted_method, ")")
  )

  p1 / p2 / p3 / p4
}


#' Write results tables + summary counts to Excel
#'
#' @param results         Wide table from add_significance_labels()
#' @param comparisons     Character vector of comparison names
#' @param output_path     Base output directory (tables/ subdir must exist)
#' @param adjusted_method Label for the summary table row
#' @return Named list: data_results, significant_adjusted_data,
#'   significant_pval_data, significant_qval_data, significant_holmval_data,
#'   and differential_list (one filtered+sorted data.frame per comparison) —
#'   the same objects data_analysis() used to assign into .GlobalEnv
export_de_tables <- function(results, comparisons, output_path, adjusted_method) {
  tables_dir <- file.path(output_path, "tables")

  significant_adjusted_data <- dplyr::filter(results, significant)
  significant_pval_data     <- dplyr::filter(results, significance_pval)
  significant_qval_data     <- dplyr::filter(results, significance_qval)
  significant_holmval_data  <- dplyr::filter(results, significance_holmval)

  writexl::write_xlsx(results, file.path(tables_dir, "data_results.xlsx"))
  writexl::write_xlsx(significant_adjusted_data, file.path(tables_dir, "significant_adjusted_data.xlsx"))
  writexl::write_xlsx(significant_pval_data, file.path(tables_dir, "significant_pval_data.xlsx"))
  writexl::write_xlsx(significant_qval_data, file.path(tables_dir, "significant_qval_data.xlsx"))
  writexl::write_xlsx(significant_holmval_data, file.path(tables_dir, "significant_holmval_data.xlsx"))

  n_significant <- data.frame(
    `N significant proteins by` = c("p-value", "adjusted p-value", "q-value", paste0(adjusted_method, " adjusted p-value")),
    `Number of proteins` = c(
      nrow(significant_pval_data), nrow(significant_adjusted_data),
      nrow(significant_qval_data), nrow(significant_holmval_data)
    ),
    check.names = FALSE
  )
  writexl::write_xlsx(n_significant, file.path(tables_dir, "number_significant_proteins_by_p-value.xlsx"))
  print(n_significant)

  differential_list <- lapply(comparisons, function(comp) {
    diff_col  <- paste0(comp, "_diffexpressed_pval")
    score_col <- paste0(comp, "_Score")
    diff <- results %>%
      dplyr::filter(.data[[diff_col]] %in% c("UP", "DOWN")) %>%
      dplyr::select(
        name, ID, genename, dplyr::all_of(diff_col),
        dplyr::all_of(paste0(comp, c("_ratio", "_p.val", "_p.adj", "_q.val", "_Score")))
      ) %>%
      dplyr::arrange(dplyr::desc(.data[[score_col]]))
    writexl::write_xlsx(diff, file.path(tables_dir, paste0(comp, "_df.xlsx")))
    diff
  })
  names(differential_list) <- comparisons

  list(
    data_results = results,
    significant_adjusted_data = significant_adjusted_data,
    significant_pval_data = significant_pval_data,
    significant_qval_data = significant_qval_data,
    significant_holmval_data = significant_holmval_data,
    differential_list = differential_list
  )
}


#' Differential expression analysis on an imputed SummarizedExperiment
#'
#' Drop-in replacement for data_analysis() (04_data_analysis.R), using limma
#' directly instead of DEP::analyze_dep(). See file header for the schema/
#' compatibility contract with 05_Plots.R and downstream scripts.
#'
#' @param se              Imputed SummarizedExperiment (assay "intensity", no NAs)
#' @param comparisons     Character vector of "A_vs_B" contrasts
#' @param output_path     Base output directory (tables/, figures/, RData/)
#' @param alpha           Significance threshold (default: global p_val)
#' @param lfc             |log2FC| threshold (default: global FC)
#' @param adjusted_method Method for the extra p.adjust() column/significance flavour
#' @return Invisibly, the list from export_de_tables(); also assigns
#'   data_results, significant_adjusted_data, significant_pval_data,
#'   significant_qval_data, significant_holmval_data, and one <comparison>_df
#'   per comparison into .GlobalEnv (matching data_analysis()'s contract)
run_de_analysis <- function(se, comparisons, output_path,
                            alpha = p_val, lfc = FC,
                            adjusted_method = c("BH", "BY", "bonferroni", "holm", "fdr", "none")) {
  adjusted_method <- match.arg(adjusted_method)
  create_directories(output_path)

  design <- build_de_design(se)
  contrasts_mx <- build_de_contrasts(design, comparisons)
  fit <- fit_de_model(se, design, contrasts_mx)

  results <- build_results_table(se, fit, comparisons, adjusted_method)
  results <- add_significance_labels(results, comparisons, alpha, lfc)

  pval_plot <- plot_pvalue_distributions(results, adjusted_method)
  save_plot("pq_value_distribution", pval_plot, subdir = "DE", width = 10, height = 12)

  save(fit, design, contrasts_mx, file = file.path(output_path, "RData", "de_fit.RData"))

  out <- export_de_tables(results, comparisons, output_path, adjusted_method)

  assign("data_results", out$data_results, envir = .GlobalEnv)
  assign("significant_adjusted_data", out$significant_adjusted_data, envir = .GlobalEnv)
  assign("significant_pval_data", out$significant_pval_data, envir = .GlobalEnv)
  assign("significant_qval_data", out$significant_qval_data, envir = .GlobalEnv)
  assign("significant_holmval_data", out$significant_holmval_data, envir = .GlobalEnv)
  for (comp in comparisons) {
    assign(paste0(comp, "_df"), out$differential_list[[comp]], envir = .GlobalEnv)
  }

  invisible(out)
}
