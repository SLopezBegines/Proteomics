# Save Plots
#' Save a ggplot to TIFF and PDF; increment global image_number counter
#' @param plotname  Base name of the file (no extension)
#' @param plot      A ggplot object
#' @param width     Width in inches
#' @param height    Height in inches
#' @param subdir    Subdirectory under figures (e.g. "QC", "DE")
save_plot <- function(plotname, plot, width = 8, height = 6, subdir = "", print = TRUE) {
  dir <- file.path(params$output_path, "figures", subdir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  filename <- paste0(dir, "/", sprintf("%03d", params$image_number), "_", plotname)
  tryCatch(
    {
      if (inherits(plot, c("ggplot", "patchwork"))) {
        ggsave(paste0(filename, params$tiff_extension), plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
        ggsave(paste0(filename, params$pdf_extension), plot, width = width, height = height, units = "in", bg = "white")
      } else if (inherits(plot, c("Heatmap", "HeatmapList", "UpSet"))) {
        # ComplexHeatmap objects must be rendered via draw() inside the device
        tiff(paste0(filename, params$tiff_extension), width = width, height = height, units = "in", res = 300)
        ComplexHeatmap::draw(plot)
        dev.off()
        pdf(paste0(filename, params$pdf_extension), width = width, height = height)
        ComplexHeatmap::draw(plot)
        dev.off()
      } else if (is.function(plot)) {
        # base-graphics functions that draw directly to the device (e.g. DimHeatmap)
        tiff(paste0(filename, params$tiff_extension), width = width, height = height, units = "in", res = 300)
        plot()
        dev.off()
        pdf(paste0(filename, params$pdf_extension), width = width, height = height)
        plot()
        dev.off()
      } else {
        warning(sprintf("[SAVE_PLOT] '%s' — unsupported plot type '%s', skipping.", plotname, class(plot)[1]))
        return(invisible(NULL))
      }
      params$image_number <<- params$image_number + 1
    },
    error = function(e) {
      warning(sprintf("[SAVE_PLOT] Failed to save '%s': %s", filename, e$message))
    }
  )
  if (print) {
    print(plot)
  }
}


# Prepare DE gene sets ####
# Reusable helper to build per-comparison gene-set lists (UP/DOWN/UPDOWN)
# from a `data_results` table. Used as input preparation for GO, STRING,
# gseGO, gseKEGG, etc.

prepare_DE_gene_sets <- function(data_results,
                                 comparisons,
                                 direction = c("UP", "DOWN"),
                                 independent_UPDOWN = TRUE,
                                 id_cols = "ID",
                                 flatten = TRUE,
                                 diffexpressed_type = c("qval", "pval", "padj_holm_pval"),
                                 top_n = NULL,
                                 keep_stats = FALSE) {
  #' @param data_results Data frame produced by data_analysis() (04_data_analysis.R),
  #'   with `<comparison>_diffexpressed_<diffexpressed_type>` columns
  #'   ("UP"/"DOWN"/"NO") for each comparison in `comparisons`.
  #' @param comparisons Character vector of comparison names (e.g. "CTRL_vs_WT").
  #' @param direction Which categories to keep: "UP", "DOWN", or both.
  #' @param independent_UPDOWN If TRUE, UP and DOWN gene sets are returned
  #'   separately (one element per comparison and direction); if FALSE, a
  #'   single combined set per comparison is returned.
  #' @param id_cols Column(s) of `data_results` to keep for each gene set.
  #' @param flatten If TRUE, each one-column gene set is unlisted into a
  #'   plain named vector (input format expected by enrichGO). If FALSE,
  #'   each gene set is kept as a data frame (e.g. for STRINGdb).
  #' @param diffexpressed_type Which significance criterion to use for the
  #'   "diffexpressed" column: "qval" (Storey's q-value, default), "pval"
  #'   (raw p-value), or "padj_holm_pval" (adjusted p-value via `adjusted_method`).
  #' @param top_n If not NULL, keep only the `top_n` proteins with the largest
  #'   absolute fold change (`<comparison>_ratio`) in each set: the highest
  #'   for UP, the most negative for DOWN, largest |FC| overall for UPDOWN.
  #' @param keep_stats If TRUE, add `foldchange` (`<comparison>_ratio`) and
  #'   `pvalue` (`<comparison>_p.val`) columns to each set (needed by gseGO /
  #'   gseKEGG). Requires `flatten = FALSE`.
  #' @return A flat named list of gene sets (vectors if `flatten = TRUE`,
  #'   data frames otherwise), named "names_<comparison>_UP"/"_DN"/"_UPDOWN".

  direction <- match.arg(direction, choices = c("UP", "DOWN"), several.ok = TRUE)
  diffexpressed_type <- match.arg(diffexpressed_type)
  if (keep_stats && flatten) stop("keep_stats = TRUE requires flatten = FALSE")

  # Optionally keep the top_n largest |FC|, then select the requested columns
  build_set <- function(sub_df, comp) {
    ratio_col <- paste0(comp, "_ratio")
    if (!is.null(top_n)) {
      sub_df <- sub_df[order(abs(sub_df[[ratio_col]]), decreasing = TRUE, na.last = NA), , drop = FALSE]
      sub_df <- utils::head(sub_df, top_n)
    }
    out <- as.data.frame(sub_df[, id_cols, drop = FALSE])
    if (keep_stats) {
      out$foldchange <- sub_df[[ratio_col]]
      out$pvalue <- sub_df[[paste0(comp, "_p.val")]]
    }
    out
  }

  gene_sets <- lapply(comparisons, function(comp) {
    diffexpressed_var <- paste0(comp, "_diffexpressed_", diffexpressed_type)

    if (independent_UPDOWN) {
      sets <- lapply(direction, function(dir) {
        sub_df <- data_results %>% dplyr::filter(!!sym(diffexpressed_var) == dir)
        build_set(sub_df, comp)
      })
      names(sets) <- paste0("names_", comp, "_", ifelse(direction == "UP", "UP", "DN"))
    } else {
      sub_df <- data_results %>% dplyr::filter(!!sym(diffexpressed_var) %in% direction)
      sets <- list(build_set(sub_df, comp))
      names(sets) <- paste0("names_", comp, "_UPDOWN")
    }
    sets
  })

  # Flatten comparisons -> flat list of gene sets
  gene_sets <- unlist(gene_sets, recursive = FALSE)

  if (flatten) {
    # Flatten one-column data frames -> plain named vectors
    gene_sets <- unlist(gene_sets, recursive = FALSE)
  }

  return(gene_sets)
}


# GSEA helpers (shared by 08_gseGO.R and 09_gseKEGG.R) ####

#' Ranked fold-change vectors for GSEA
#' @param sets  Output of prepare_DE_gene_sets(..., keep_stats = TRUE, flatten = FALSE)
#' @return Named list of numeric vectors (names = `ID`), NA-free, decreasing.
prepare_gse_gene_lists <- function(sets) {
  lapply(sets, function(df) sort(stats::na.omit(stats::setNames(df$foldchange, df$ID)), decreasing = TRUE))
}

#' Lollipop plot of the top terms of a gseaResult data frame
#' @param df     `@result` of a gseaResult
#' @param title  Plot title
#' @param n_terms Number of terms to show (highest enrichmentScore)
#' @param x_label Axis label for the terms
#' @return A ggplot, or NULL if `df` has no rows. Facets by ONTOLOGY when present.
plot_gse_lolliplot <- function(df, title, n_terms = 15, x_label = "Term") {
  if (nrow(df) == 0) {
    return(NULL)
  }
  p <- df %>%
    dplyr::mutate(Description = stats::reorder(Description, enrichmentScore)) %>%
    dplyr::slice_max(enrichmentScore, n = n_terms, with_ties = FALSE) %>%
    ggplot(aes(x = Description, y = enrichmentScore, colour = p.adjust)) +
    geom_segment(aes(xend = Description, y = 0, yend = enrichmentScore)) +
    geom_point(aes(size = setSize - 1), show.legend = TRUE) +
    scale_color_viridis_c(option = "viridis", direction = 1) +
    coord_flip() +
    theme_minimal() +
    labs(size = "N. of genes", x = x_label, y = "enrichmentScore") +
    ggtitle(title)
  if ("ONTOLOGY" %in% names(df)) p <- p + facet_wrap(ONTOLOGY ~ ., scales = "free_y", ncol = 1)
  p
}
