# de_plots.R
#
# DEP-free visualization functions for differential-expression results.
#
# Ports every plot from 05_Plots.R into standalone functions, replacing all
# DEP:: calls (DEP::plot_cor, DEP::plot_heatmap, DEP::plot_volcano,
# DEP::plot_single, DEP::theme_DEP1) — the DEP-based PCA plot and
# DEP::plot_volcano were dropped rather than ported: PCA is already covered
# by qc_proteomics.R::plot_pca()/plot_pca_comparison(), and the ggplot
# q-value volcano covers the same need as DEP::plot_volcano() without DEP.
# All functions work directly off the wide `data_results` table produced by
# run_de_analysis() (code/de_analysis.R) — lfq_* intensity columns,
# <comparison>_ratio/_q.val/_diffexpressed_qval/_dif_label_qval, `name`,
# `significant` — so nothing here needs to reload a model object.
#
# Function index:
#   .too_few_rows()          [internal]
#   .safe_cutree_k()         [internal]
#   plot_sample_correlation()
#   plot_significant_heatmap()
#   plot_protein_barplot()
#   plot_protein_barplots_batched()
#   plot_pca_loadings()
#   plot_volcano_qval()
#   plot_fc_correlation()
#   plot_fc_correlations_all()
#   plot_ratio_heatmap()
#   plot_ratio_heatmap_excluded()
#   plot_intensity_heatmap()
#   plot_intensity_heatmap_excluded()


#' Warn and signal "skip" when there isn't enough data to cluster/correlate
#'
#' Several plots below filter data_results down to `significant`/
#' `significance_qval` proteins first; when that filter leaves too few rows,
#' stats::cor()/stats::hclust() error out with cryptic messages ("'x' is
#' empty", "must have n >= 2 objects to cluster"). This turns that into an
#' informative warning and a clean early return instead.
#'
#' @param n        Number of rows actually available
#' @param min_n    Minimum required
#' @param fun_name Calling function's name, for the warning message
#' @return TRUE if the caller should skip plotting (not enough rows), FALSE otherwise
.too_few_rows <- function(n, min_n, fun_name) {
  if (n < min_n) {
    warning(sprintf(
      "[%s] Only %d protein(s) available (need >= %d) — likely too few significant proteins at the current thresholds. Skipping.",
      fun_name, n, min_n
    ), call. = FALSE)
    return(TRUE)
  }
  FALSE
}


#' Cap a requested cluster count at what stats::cutree() can actually do (1..n-1)
#' @param k Requested number of clusters
#' @param n Number of rows being clustered
#' @return Integer between 1 and n-1
.safe_cutree_k <- function(k, n) {
  max(1L, min(k, n - 1L))
}


#' Pearson correlation heatmap across samples
#'
#' Replaces DEP::plot_cor(dep_analysis, significant = TRUE, lower = 0, upper = 1, pal = "Reds").
#'
#' @param data_results     Wide table with lfq_* columns (run_de_analysis() output)
#' @param significant_only Restrict the correlation to `significant` proteins
#'   (default TRUE, matches DEP's default)
#' @return ComplexHeatmap::Heatmap, or NULL (with a warning) if fewer than 2
#'   proteins pass the filter
plot_sample_correlation <- function(data_results, significant_only = TRUE) {
  mat <- data_results
  if (significant_only) mat <- dplyr::filter(mat, significant)

  lfq_cols <- grep("^lfq_", names(mat), value = TRUE)
  intensity <- as.matrix(mat[lfq_cols])
  colnames(intensity) <- sub("^lfq_", "", colnames(intensity))

  if (.too_few_rows(nrow(intensity), 2, "plot_sample_correlation")) return(invisible(NULL))

  cor_mat <- stats::cor(intensity, use = "pairwise.complete.obs")

  ComplexHeatmap::Heatmap(
    cor_mat,
    name = "Pearson r",
    col = circlize::colorRamp2(range(cor_mat), c("white", "#B2182B")),
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid::grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = grid::gpar(fontsize = 8))
    },
    column_title = "Sample correlation"
  )
}


#' Row-centered heatmap of significant proteins, k-means clustered
#'
#' Replaces DEP::plot_heatmap(dep_analysis, type = "centered", kmeans = TRUE,
#' k = 2, col_limit = 2, indicate = c("condition", "replicate")).
#'
#' @param data_results Wide table with lfq_* columns and a `significant` flag
#' @param exp_design   Experimental design data.frame (columns_to_rename, condition, replicate)
#' @param k            Number of k-means clusters (row_km)
#' @param col_limit    Clip centered values to [-col_limit, col_limit]
#' @return ComplexHeatmap::Heatmap, or NULL (with a warning) if fewer than 2
#'   significant proteins are found (row clustering needs >= 2)
plot_significant_heatmap <- function(data_results, exp_design, k = 2, col_limit = 2) {
  mat <- dplyr::filter(data_results, significant)

  lfq_cols <- grep("^lfq_", names(mat), value = TRUE)
  intensity <- as.matrix(mat[lfq_cols])
  rownames(intensity) <- mat$name
  colnames(intensity) <- sub("^lfq_", "", colnames(intensity))

  if (.too_few_rows(nrow(intensity), 2, "plot_significant_heatmap")) return(invisible(NULL))

  centered <- t(scale(t(intensity), center = TRUE, scale = FALSE))
  centered[centered >  col_limit] <-  col_limit
  centered[centered < -col_limit] <- -col_limit

  sample_order <- match(colnames(centered), exp_design$columns_to_rename)
  col_anno <- ComplexHeatmap::HeatmapAnnotation(
    condition = exp_design$condition[sample_order],
    replicate = as.character(exp_design$replicate[sample_order])
  )

  # k-means needs at least k points; fall back to plain hierarchical
  # clustering (no split) when fewer significant proteins than k
  row_km <- if (nrow(centered) >= k) k else NULL

  ComplexHeatmap::Heatmap(
    centered,
    name = "Centered\nlog2 intensity",
    col = circlize::colorRamp2(c(-col_limit, 0, col_limit), c("#2166AC", "white", "#B2182B")),
    top_annotation = col_anno,
    row_km = row_km,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 6),
    column_title = "Significant proteins (row-centered)"
  )
}


#' Grouped intensity barplot for a set of proteins, one panel per protein
#'
#' Replaces DEP::plot_single(dep_analysis, proteins = ..., type = "centered").
#'
#' @param data_results Wide table with `name` and lfq_* columns
#' @param exp_design   Experimental design data.frame (columns_to_rename, condition)
#' @param proteins     Character vector of protein names (data_results$name) to plot
#' @return ggplot, facetted by protein
plot_protein_barplot <- function(data_results, exp_design, proteins) {
  lfq_cols <- grep("^lfq_", names(data_results), value = TRUE)

  df <- data_results %>%
    dplyr::filter(name %in% proteins) %>%
    dplyr::select(name, dplyr::all_of(lfq_cols)) %>%
    tidyr::pivot_longer(-name, names_to = "sample", values_to = "intensity") %>%
    dplyr::mutate(sample = sub("^lfq_", "", sample))

  df$condition <- exp_design$condition[match(df$sample, exp_design$columns_to_rename)]

  ggplot2::ggplot(df, ggplot2::aes(x = condition, y = intensity, fill = condition)) +
    ggplot2::stat_summary(fun = mean, geom = "bar") +
    ggplot2::geom_jitter(width = 0.15, size = 1) +
    ggplot2::facet_wrap(~name, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = "log2 intensity", x = NULL) +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' Batch a protein list into pages and barplot each page
#'
#' Adapted from 05_Plots.R's "Barplots" loop, which called DEP::plot_single()
#' `per_page` proteins at a time. Unlike the original loop (which used
#' integer division `length(proteins) / per_page` and silently dropped a
#' partial last page), this uses ceiling() so every protein is included.
#'
#' @param data_results Wide table with `name` and lfq_* columns
#' @param exp_design   Experimental design data.frame
#' @param proteins     Character vector of protein names to plot
#' @param per_page     Proteins per page (default 16, matches the original)
#' @return List of ggplots, one per page
plot_protein_barplots_batched <- function(data_results, exp_design, proteins, per_page = 16) {
  n_pages <- ceiling(length(proteins) / per_page)
  lapply(seq_len(n_pages), function(page) {
    start <- (page - 1) * per_page + 1
    end   <- min(page * per_page, length(proteins))
    plot_protein_barplot(data_results, exp_design, proteins[start:end])
  })
}


#' PCA biplot of proteins, coloured by a significance column
#'
#' Adapted from 05_Plots.R's "PCA plot loadings" block. PCA is fit on
#' significant/significance_qval + all lfq_* columns, one row per protein
#' (this mirrors the original exactly; note it PCAs proteins using their
#' significance flags as extra "features" alongside intensities, not a
#' clean sample-space PCA — kept as-is since that's what the legacy plot did).
#'
#' The original also set `pca_res$repetition <- rep(1:4, times = 4)`
#' (hardcoded to 16 samples) after the prcomp() call; it was never read by
#' autoplot() and doesn't match this dataset's actual sample count (26), so
#' it's dropped here as dead code.
#'
#' @param data_results Wide table with `name`, `significant`,
#'   `significance_qval`, and lfq_* columns
#' @param colour_by    Column of `data_results` to colour points by
#'   (e.g. "significant" or "significance_qval")
#' @return ggplot (via ggfortify::autoplot.prcomp)
plot_pca_loadings <- function(data_results, colour_by = "significant") {
  lfq_cols <- grep("^lfq_", names(data_results), value = TRUE)
  pca_matrix <- as.matrix(data_results[, c("significant", "significance_qval", lfq_cols)])
  rownames(pca_matrix) <- data_results$name

  pca_res <- stats::prcomp(pca_matrix, scale. = FALSE)

  ggplot2::autoplot(pca_res,
    data = data_results, colour = colour_by, shape = TRUE,
    label = FALSE, label.size = 3,
    loadings = FALSE, loadings.colour = "blue",
    loadings.label = TRUE, loadings.label.size = 3,
    frame = FALSE, frame.type = "norm"
  )
}


#' ggplot theme used by plot_volcano_qval()
#'
#' Reimplements DEP::theme_DEP1() (a plain ggplot2::theme_bw() variant) so
#' this file has no DEP dependency at all, not even a cosmetic one.
#' @return ggplot2 theme
theme_volcano <- function() {
  basesize <- 12
  theme <- ggplot2::theme_bw(base_size = basesize)
  theme$plot.title$face <- "bold"
  theme$plot.title$size <- basesize + 2
  theme$plot.title$hjust <- 0.5
  theme$axis.title.x$size <- basesize + 2
  theme$axis.title.y$size <- basesize + 2
  theme$axis.text$size <- basesize
  theme$axis.text$colour <- "black"
  theme$legend.title$size <- basesize + 2
  theme$legend.text$size <- basesize
  theme$strip.text$face <- "bold"
  theme$strip.text$size <- basesize + 2
  theme$strip.text$colour <- "black"
  theme
}


#' Volcano plot (log2FC vs -log10(q-value)) for one comparison
#'
#' Replaces DEP::plot_volcano(dep_analysis, contrast = ..., adjusted = FALSE)
#' — adapted from 05_Plots.R's "ggplot q-value" volcano block, which already
#' covered the same need without DEP; this just makes it a function.
#'
#' @param data_results Wide table with <comparison>_ratio/_q.val/
#'   _diffexpressed_qval/_dif_label_qval columns
#' @param comparison   Single comparison name in "A_vs_B" form (e.g.
#'   "CTRL_vs_WT"); the contrast is A - B, so A is annotated on the right
#'   (log2FC > 0, enriched in A) and B on the left (log2FC < 0, enriched in B)
#' @param p_val        Horizontal significance threshold line (q-value)
#' @param fc           Vertical |log2FC| threshold lines
#' @return ggplot
plot_volcano_qval <- function(data_results, comparison, p_val, fc) {
  ratio <- paste0(comparison, "_ratio")
  qval  <- paste0(comparison, "_q.val")
  diff_ <- paste0(comparison, "_diffexpressed_qval")
  label <- paste0(comparison, "_dif_label_qval")

  mycolors <- c(DOWN = "#F39B7F", UP = "#4DBBD5", NO = "#A6A6A6")
  groups <- strsplit(comparison, "_vs_", fixed = TRUE)[[1]]

  p <- suppressWarnings(
    ggplot2::ggplot(data_results, ggplot2::aes(
      x = .data[[ratio]], y = -log10(.data[[qval]]), col = .data[[diff_]]
    )) +
      ggplot2::geom_point() +
      theme_volcano() +
      ggplot2::ggtitle(paste("Volcano", comparison)) +
      ggplot2::labs(
        x = paste("Fold-Change", comparison),
        y = expression(paste("-log"[10], "(q-value)"))
      ) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_vline(xintercept = c(-fc, fc), col = "red", linetype = "dashed") +
      ggplot2::geom_hline(yintercept = -log10(p_val), col = "red", linetype = "dashed") +
      ggplot2::scale_colour_manual(values = mycolors) +
      ggrepel::geom_text_repel(ggplot2::aes(label = .data[[label]]), max.overlaps = 15) +
      ggplot2::guides(colour = "none")
  )

  # Direction labels: A (right, log2FC > 0) vs B (left, log2FC < 0)
  if (length(groups) == 2) {
    p <- p +
      ggplot2::annotate("text",
        x = Inf, y = Inf, label = groups[1], hjust = 1.1, vjust = 1.5,
        fontface = "bold", size = 4
      ) +
      ggplot2::annotate("text",
        x = -Inf, y = Inf, label = groups[2], hjust = -0.1, vjust = 1.5,
        fontface = "bold", size = 4
      )
  }

  p
}


#' Scatter plot comparing log2FC between two comparisons
#'
#' Adapted from 05_Plots.R's "Fold-Change correlations" block.
#'
#' @param data_results  Wide table with <comparison>_ratio columns
#' @param comparison_x  First comparison name
#' @param comparison_y  Second comparison name
#' @return ggplot
plot_fc_correlation <- function(data_results, comparison_x, comparison_y) {
  x_col <- paste0(comparison_x, "_ratio")
  y_col <- paste0(comparison_y, "_ratio")

  ggplot2::ggplot(data_results, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]])) +
    ggplot2::geom_point(shape = 21, fill = "white", size = 2) +
    ggrepel::geom_text_repel(ggplot2::aes(label = name), max.overlaps = 15) +
    ggplot2::labs(title = paste(x_col, "vs", y_col), x = x_col, y = y_col)
}


#' All pairwise log2FC correlation plots across comparisons
#'
#' @param data_results Wide table with <comparison>_ratio columns
#' @param comparisons  Character vector of >= 2 comparison names
#' @return Named list of ggplots, one per pair ("compA__vs__compB"); empty
#'   list (with a warning) if fewer than 2 comparisons are given
plot_fc_correlations_all <- function(data_results, comparisons) {
  if (length(comparisons) < 2) {
    warning("Need at least 2 comparisons for FC correlation plots.")
    return(list())
  }
  pairs <- utils::combn(comparisons, 2, simplify = FALSE)
  plots <- lapply(pairs, function(p) plot_fc_correlation(data_results, p[1], p[2]))
  names(plots) <- vapply(pairs, function(p) paste(p[1], p[2], sep = "__vs__"), character(1))
  plots
}


#' Heatmap of log2FC ratios for filtered proteins, one column per comparison
#'
#' Adapted from 05_Plots.R's "HeatMap padj"/"HeatMap pval" blocks (both used
#' the same pheatmap logic, just filtered proteins by a different
#' significance column — generalised here into one function).
#'
#' @param data_results Wide table with <comparison>_ratio columns
#' @param comparisons  Character vector of comparison names (also sets column order)
#' @param filter_col   Logical column to filter proteins by (e.g. "significant", "significance_qval")
#' @param cutree_rows  Number of row clusters to cut the dendrogram into
#' @param title        pheatmap main title
#' @return pheatmap object, or NULL (with a warning) if no protein passes `filter_col`
plot_ratio_heatmap <- function(data_results, comparisons, filter_col = "significant",
                                cutree_rows = 4, title = "Fold-change ratio") {
  ratio_cols <- paste0(comparisons, "_ratio")

  df <- data_results[data_results[[filter_col]] == TRUE, c("name", ratio_cols)]
  mat <- as.matrix(df[, ratio_cols, drop = FALSE])
  rownames(mat) <- df$name

  if (.too_few_rows(nrow(mat), 1, "plot_ratio_heatmap")) return(invisible(NULL))

  # hclust/cutree need >= 2 rows; fall back to unclustered for a single row
  can_cluster <- nrow(mat) >= 2

  pheatmap::pheatmap(mat,
    cutree_rows = if (can_cluster) .safe_cutree_k(cutree_rows, nrow(mat)) else 1,
    cluster_cols = FALSE,
    cluster_rows = can_cluster,
    scale = "none",
    annotation_legend = TRUE,
    labels_col = comparisons,
    labels_row = df$name,
    main = title,
    display_numbers = TRUE
  )
}


#' Ratio heatmap, hierarchically clustered, with selected clusters dropped
#'
#' Adapted from 05_Plots.R's "HeatMap pval excluded" block. Which clusters to
#' drop is a visual/manual call made after inspecting the dendrogram — the
#' legacy script hardcoded "clusters 1 and 5 have many uninformative
#' proteins" for its dataset. Re-tune `exclude_clusters` per dataset/run;
#' the default (none excluded) is intentionally conservative.
#'
#' @param data_results     Wide table with <comparison>_ratio columns
#' @param comparisons      Character vector of comparison names (also sets column order)
#' @param filter_col       Logical column to filter proteins by
#' @param k                Number of hclust clusters to cut into
#' @param exclude_clusters Cluster numbers (from cutree) to drop before plotting
#' @param title            pheatmap main title
#' @return list(plot = pheatmap object, kept_data = data_results rows for the
#'   surviving, non-excluded proteins, for the caller to export/inspect) —
#'   or NULL (with a warning) if there's too little data to cluster, either
#'   before or after exclusion
plot_ratio_heatmap_excluded <- function(data_results, comparisons, filter_col = "significance_qval",
                                         k = 6, exclude_clusters = integer(0),
                                         title = "Fold-change ratio") {
  ratio_cols <- paste0(comparisons, "_ratio")

  df <- data_results[data_results[[filter_col]] == TRUE, c("name", ratio_cols)]
  mat <- as.matrix(df[, ratio_cols, drop = FALSE])
  rownames(mat) <- df$name

  if (.too_few_rows(nrow(mat), 2, "plot_ratio_heatmap_excluded")) return(invisible(NULL))

  k <- .safe_cutree_k(k, nrow(mat))
  hc <- stats::hclust(stats::dist(mat))
  clusters <- stats::cutree(hc, k = k)
  mat <- mat[!(clusters %in% exclude_clusters), , drop = FALSE]

  if (.too_few_rows(nrow(mat), 1, "plot_ratio_heatmap_excluded (after exclusion)")) return(invisible(NULL))
  can_cluster <- nrow(mat) >= 2

  p <- pheatmap::pheatmap(mat,
    cutree_rows = if (can_cluster) .safe_cutree_k(k, nrow(mat)) else 1,
    cluster_cols = FALSE,
    cluster_rows = can_cluster,
    scale = "none",
    annotation_legend = TRUE,
    labels_col = comparisons,
    labels_row = rownames(mat),
    main = title,
    display_numbers = TRUE
  )

  list(
    plot = p,
    kept_data = data_results[data_results$name %in% rownames(mat), ]
  )
}


#' Heatmap of log2 LFQ intensities for filtered proteins, all samples
#'
#' Adapted from 05_Plots.R's "Log2 Intensity"/"Log2 Intensity significant"
#' blocks. The original passed a non-existent `row_names` argument to
#' pheatmap() (silently swallowed by its `...`, so row names were always
#' shown regardless of the intended TRUE/FALSE) — fixed here to the real
#' parameter, `show_rownames`.
#'
#' @param data_results   Wide table with lfq_* columns
#' @param filter_col     Logical column to filter proteins by
#' @param sample_labels  Column labels for samples (defaults to lfq_ column
#'   names with the "lfq_" prefix stripped)
#' @param cutree_rows    Number of row clusters to cut the dendrogram into
#' @param show_rownames  Show protein names on rows
#' @return pheatmap object, or NULL (with a warning) if no protein passes `filter_col`
plot_intensity_heatmap <- function(data_results, filter_col = "significant",
                                    sample_labels = NULL, cutree_rows = 4,
                                    show_rownames = TRUE) {
  lfq_cols <- grep("^lfq_", names(data_results), value = TRUE)
  if (is.null(sample_labels)) sample_labels <- sub("^lfq_", "", lfq_cols)

  df <- data_results[data_results[[filter_col]] == TRUE, c("name", lfq_cols)]
  mat <- as.matrix(df[, lfq_cols, drop = FALSE])
  rownames(mat) <- df$name

  if (.too_few_rows(nrow(mat), 1, "plot_intensity_heatmap")) return(invisible(NULL))

  # hclust/cutree need >= 2 rows; fall back to unclustered for a single row
  can_cluster <- nrow(mat) >= 2

  pheatmap::pheatmap(mat,
    cutree_rows = if (can_cluster) .safe_cutree_k(cutree_rows, nrow(mat)) else 1,
    cluster_cols = FALSE,
    cluster_rows = can_cluster,
    scale = "none",
    show_rownames = show_rownames,
    annotation_legend = TRUE,
    labels_col = sample_labels,
    main = "Log2 Intensity"
  )
}


#' Intensity heatmap, hierarchically clustered, with selected clusters dropped
#'
#' Adapted from 05_Plots.R's "Log2 Intensity Excluded" block. Same manual
#' cluster-exclusion caveat as plot_ratio_heatmap_excluded() — re-tune
#' `exclude_clusters` per dataset/run after inspecting the dendrogram.
#'
#' @param data_results     Wide table with lfq_* columns
#' @param filter_col       Logical column to filter proteins by
#' @param sample_labels    Column labels for samples
#' @param k                Number of hclust clusters to cut into
#' @param exclude_clusters Cluster numbers (from cutree) to drop before plotting
#' @return list(plot = pheatmap object, kept_data = data_results rows for the
#'   surviving, non-excluded proteins) — or NULL (with a warning) if there's
#'   too little data to cluster, either before or after exclusion
plot_intensity_heatmap_excluded <- function(data_results, filter_col = "significance_qval",
                                             sample_labels = NULL, k = 6,
                                             exclude_clusters = integer(0)) {
  lfq_cols <- grep("^lfq_", names(data_results), value = TRUE)
  if (is.null(sample_labels)) sample_labels <- sub("^lfq_", "", lfq_cols)

  df <- data_results[data_results[[filter_col]] == TRUE, c("name", lfq_cols)]
  mat <- as.matrix(df[, lfq_cols, drop = FALSE])
  rownames(mat) <- df$name

  if (.too_few_rows(nrow(mat), 2, "plot_intensity_heatmap_excluded")) return(invisible(NULL))

  k <- .safe_cutree_k(k, nrow(mat))
  hc <- stats::hclust(stats::dist(mat))
  clusters <- stats::cutree(hc, k = k)
  mat <- mat[!(clusters %in% exclude_clusters), , drop = FALSE]

  if (.too_few_rows(nrow(mat), 1, "plot_intensity_heatmap_excluded (after exclusion)")) return(invisible(NULL))
  can_cluster <- nrow(mat) >= 2

  p <- pheatmap::pheatmap(mat,
    cutree_rows = if (can_cluster) .safe_cutree_k(k, nrow(mat)) else 1,
    cluster_cols = FALSE,
    cluster_rows = can_cluster,
    scale = "none",
    show_rownames = TRUE,
    annotation_legend = TRUE,
    labels_col = sample_labels,
    main = "Log2 Intensity"
  )

  list(
    plot = p,
    kept_data = data_results[data_results$name %in% rownames(mat), ]
  )
}
