# Save Plots
#' Save a ggplot to TIFF and PDF; increment global image_number counter
#' @param plotname  Base name of the file (no extension)
#' @param plot      A ggplot object
#' @param width     Width in inches
#' @param height    Height in inches
#' @param subdir    Subdirectory under figures (e.g. "QC", "DE")
save_plot <- function(plotname, plot, width = 8, height = 6, subdir = "", print = TRUE) {
  dir <- file.path(output_path, "figures", subdir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  filename <- paste0(dir, "/", sprintf("%03d", image_number), "_", plotname)
  tryCatch(
    {
      if (inherits(plot, c("ggplot", "patchwork"))) {
        ggsave(paste0(filename, tiff_extension), plot, width = width, height = height, units = "in", dpi = 300, bg = "white")
        ggsave(paste0(filename, pdf_extension), plot, width = width, height = height, units = "in", bg = "white")
      } else if (inherits(plot, c("Heatmap", "HeatmapList", "UpSet"))) {
        # ComplexHeatmap objects must be rendered via draw() inside the device
        tiff(paste0(filename, tiff_extension), width = width, height = height, units = "in", res = 300)
        ComplexHeatmap::draw(plot)
        dev.off()
        pdf(paste0(filename, pdf_extension), width = width, height = height)
        ComplexHeatmap::draw(plot)
        dev.off()
      } else if (is.function(plot)) {
        # base-graphics functions that draw directly to the device (e.g. DimHeatmap)
        tiff(paste0(filename, tiff_extension), width = width, height = height, units = "in", res = 300)
        plot()
        dev.off()
        pdf(paste0(filename, pdf_extension), width = width, height = height)
        plot()
        dev.off()
      } else {
        warning(sprintf("[SAVE_PLOT] '%s' — unsupported plot type '%s', skipping.", plotname, class(plot)[1]))
        return(invisible(NULL))
      }
      image_number <<- image_number + 1
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
                                 diffexpressed_type = c("qval", "pval", "padj_holm_pval")) {
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
  #' @return A flat named list of gene sets (vectors if `flatten = TRUE`,
  #'   data frames otherwise), named "names_<comparison>_UP"/"_DN"/"_UPDOWN".

  direction <- match.arg(direction, choices = c("UP", "DOWN"), several.ok = TRUE)
  diffexpressed_type <- match.arg(diffexpressed_type)

  gene_sets <- lapply(comparisons, function(comp) {
    diffexpressed_var <- paste0(comp, "_diffexpressed_", diffexpressed_type)

    if (independent_UPDOWN) {
      sets <- lapply(direction, function(dir) {
        sub_df <- data_results %>% dplyr::filter(!!sym(diffexpressed_var) == dir)
        as.data.frame(sub_df[, id_cols, drop = FALSE])
      })
      names(sets) <- paste0("names_", comp, "_", ifelse(direction == "UP", "UP", "DN"))
    } else {
      sub_df <- data_results %>% dplyr::filter(!!sym(diffexpressed_var) %in% direction)
      sets <- list(as.data.frame(sub_df[, id_cols, drop = FALSE]))
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
