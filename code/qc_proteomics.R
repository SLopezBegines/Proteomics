# qc_proteomics.R
#
# Modular proteomics QC, normalization, and imputation functions.
# All functions are pure: they return their result and produce no side effects
# (no writing to disk, no global assignments). Saving is the caller's responsibility.
#
# Does NOT depend on DEP. Uses:
#   SummarizedExperiment, S4Vectors  (Bioc — core infra)
#   PRONE                            (Bioc — normalization; install: renv::install("daisybio/PRONE"))
#   ComplexHeatmap, RColorBrewer     (Bioc — missing-value heatmap)
#   ggplot2, dplyr, tidyr, purrr, broom, patchwork, ggrepel, scales
#   truncnorm                        (CRAN — MNAR imputation)
#   VIM                              (CRAN — kNN MAR imputation)
#   rrcov                            (CRAN — robust PCA outlier detection)
#
# Function index:
#   1. SE CONSTRUCTION
#      create_se()
#   2. QC — PRE-FILTERING
#      plot_id_overlap()
#      plot_upset()
#      filter_missing()
#      filter_missing_frac()
#      plot_sample_numbers()
#      plot_protein_coverage()
#      plot_missing_heatmap()
#      plot_intensity_detect()
#   3. NORMALIZATION — PRONE
#      normalize_prone()
#      plot_normalization_prone()
#      select_normalization()
#   4. QC — POST-NORMALIZATION
#      plot_pca()
#      plot_mean_sd()
#   5. MNAR / MAR CLASSIFICATION
#      classify_mnar_mar()
#   6. IMPUTATION — CUSTOM R (no MSnbase)
#      impute_knn_mar()
#      impute_mnar_truncnorm()
#      impute_mixed_split()
#   7. IMPUTATION QC
#      plot_imputation_distributions()
#      plot_sd_before_after()
#      plot_mnar_distribution()
#      .compute_sd_per_protein_condition()   [internal]
#   8. PCA — OUTLIER DETECTION
#      run_pca_hubert()
#      plot_pca_comparison()


# ── 1. SE CONSTRUCTION ────────────────────────────────────────────────────────

#' Build a SummarizedExperiment from a proteomics data.frame
#'
#' Replaces DEP::make_unique() + DEP::make_se().
#' Row names are taken from `genename`; duplicates are resolved by appending the
#' protein accession (`ID`). Zero values (MaxQuant missing-detection sentinel)
#' are converted to NA before log2 transformation.
#'
#' @param prot_data  data.frame with columns `genename`, `ID`, and intensity
#'                   columns whose names match `Exp_design$label`
#' @param Exp_design data.frame with columns `label`, `condition`, `replicate`
#' @return SummarizedExperiment with assay "intensity" (log2 scale)
create_se <- function(prot_data, Exp_design) {
  # Build unique row names: prefer genename, fall back to ID
  row_names <- ifelse(
    is.na(prot_data$genename) | prot_data$genename == "",
    prot_data$ID,
    prot_data$genename
  )
  # Resolve duplicates by appending the protein accession
  dup_idx <- duplicated(row_names) | duplicated(row_names, fromLast = TRUE)
  row_names[dup_idx] <- paste0(row_names[dup_idx], "_", prot_data$ID[dup_idx])
  # Final safety pass with make.unique
  row_names <- make.unique(row_names, sep = "_")

  # Extract intensity matrix in the column order defined by Exp_design
  intensity_cols <- match(Exp_design$label, colnames(prot_data))
  if (any(is.na(intensity_cols))) {
    stop(
      "Columns not found in prot_data: ",
      paste(Exp_design$label[is.na(intensity_cols)], collapse = ", ")
    )
  }
  mat <- as.matrix(prot_data[, intensity_cols])
  rownames(mat) <- row_names
  colnames(mat) <- Exp_design$columns_to_rename

  # Replace MaxQuant's zero sentinel with NA before log2 transformation
  mat[mat == 0] <- NA


  # Log2-transform only if values look like raw intensities (median > 100)
  # if (median(mat, na.rm = TRUE) > 100) {
  #  mat                 <- log2(mat)
  #  mat[is.infinite(mat)] <- NA  # guard against any residual zeros
  # }


  # Build colData from experimental design
  col_data <- S4Vectors::DataFrame(
    label     = Exp_design$columns_to_rename,
    condition = Exp_design$condition,
    replicate = Exp_design$replicate,
    row.names = Exp_design$columns_to_rename
  )

  # Carry protein-level metadata as rowData
  row_data <- S4Vectors::DataFrame(
    genename = prot_data$genename,
    ID = prot_data$ID,
    row.names = row_names
  )

  SummarizedExperiment::SummarizedExperiment(
    assays  = list(intensity = mat),
    colData = col_data,
    rowData = row_data
  )
}


# ── 2. QC — PRE-FILTERING ─────────────────────────────────────────────────────

#' Barplot: how many samples each protein is detected in
#'
#' Replaces DEP::plot_frequency(). Shows the distribution of protein
#' identification completeness across samples.
#'
#' @param se SummarizedExperiment
#' @return ggplot
plot_id_overlap <- function(se) {
  mat <- SummarizedExperiment::assay(se, "intensity")
  n_det <- rowSums(!is.na(mat))
  freq <- table(n_det)

  df <- data.frame(
    n_samples  = as.integer(names(freq)),
    n_proteins = as.integer(freq)
  )

  ggplot2::ggplot(df, ggplot2::aes(x = factor(n_samples), y = n_proteins)) +
    ggplot2::geom_col(fill = "steelblue", color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = n_proteins), vjust = -0.4, size = 3.5) +
    ggplot2::labs(
      title = "Protein identification overlap",
      x     = "Number of samples with detection",
      y     = "Number of proteins"
    ) +
    ggplot2::theme_minimal()
}


#' Every group and between-group EXCLUSIVE intersection, as one row per combination
#'
#' For each group alone (n = 1) and every combination of groups (n >= 2), finds the
#' proteins whose detection pattern across `comb_mat` matches that combination
#' EXACTLY (1 in the selected groups, 0 everywhere else). This is the same
#' exclusive-combination definition ComplexHeatmap::make_comb_mat() uses, so the
#' counts match the UpSet plot bars exactly — unlike Reduce(intersect, ...), which
#' would also count proteins shared with groups outside the combination.
#'
#' @param comb_mat Binary matrix (proteins x groups), as built in plot_upset_custom()
#' @return data.frame with columns Intersection, Proteins, Num_Proteins
.calculate_all_intersections <- function(comb_mat) {
  conditions <- colnames(comb_mat)
  n_groups <- length(conditions)
  row_codes <- do.call(paste0, lapply(seq_len(n_groups), function(i) comb_mat[, i]))

  intersection_data <- list()
  for (n in seq_len(n_groups)) {
    combinations <- combn(seq_len(n_groups), n, simplify = FALSE)
    intersections_n <- lapply(combinations, function(comb) {
      pattern <- rep(0L, n_groups)
      pattern[comb] <- 1L
      matches <- row_codes == paste0(pattern, collapse = "")
      proteins <- rownames(comb_mat)[matches]
      list(
        name         = paste(conditions[comb], collapse = "_AND_"),
        proteins     = paste(proteins, collapse = ", "),
        num_proteins = length(proteins)
      )
    })
    intersection_data <- c(intersection_data, intersections_n)
  }

  data.frame(
    Intersection = sapply(intersection_data, `[[`, "name"),
    Proteins = sapply(intersection_data, `[[`, "proteins"),
    Num_Proteins = sapply(intersection_data, `[[`, "num_proteins"),
    stringsAsFactors = FALSE
  )
}


#' UpSet plot of protein detection overlap across conditions
#'
#' Shows which proteins are detected in each condition and the size of every
#' intersection. Requires ComplexHeatmap (already a project dependency).
#'
#' @param se       SummarizedExperiment
#' @param min_size Minimum intersection size to display (default 5)
#' @return list:
#'   $plot  ComplexHeatmap UpSet object (use draw() or print() to render)
#'   $table data.frame with one row per group and per between-group intersection —
#'          columns Intersection, Proteins (comma-separated), Num_Proteins
plot_upset_custom <- function(se, min_size = 5) {
  mat <- SummarizedExperiment::assay(se, "intensity")
  conditions <- unique(SummarizedExperiment::colData(se)$condition)

  # Binary matrix: 1 = detected in >= 1 replicate of that condition
  comb_mat <- vapply(conditions, function(cond) {
    idx <- which(SummarizedExperiment::colData(se)$condition == cond)
    as.integer(rowSums(!is.na(mat[, idx, drop = FALSE])) > 0)
  }, integer(nrow(mat)))
  rownames(comb_mat) <- rownames(mat)
  colnames(comb_mat) <- conditions

  m <- ComplexHeatmap::make_comb_mat(comb_mat)
  m <- m[ComplexHeatmap::comb_size(m) >= min_size]

  p <- ComplexHeatmap::UpSet(
    m,
    comb_order       = order(ComplexHeatmap::comb_size(m), decreasing = TRUE),
    top_annotation   = ComplexHeatmap::upset_top_annotation(m, add_numbers = TRUE),
    right_annotation = ComplexHeatmap::upset_right_annotation(m, add_numbers = TRUE)
  )

  table <- .calculate_all_intersections(comb_mat)

  list(plot = p, table = table)
}


#' Filter proteins with too many missing values per condition
#'
#' Replaces DEP::filter_missval(). A protein is retained if it has at most
#' `thr` missing replicates in at least one condition.
#'
#' @param se  SummarizedExperiment with colData$condition
#' @param thr Maximum missing values allowed in the best-covered condition (default 1)
#' @return Filtered SummarizedExperiment
filter_missing <- function(se, thr = 1) {
  mat <- SummarizedExperiment::assay(se, "intensity")
  conditions <- unique(SummarizedExperiment::colData(se)$condition)

  pass <- apply(mat, 1, function(row) {
    # Keep protein if at least one condition has <= thr NAs
    any(vapply(conditions, function(cond) {
      idx <- which(SummarizedExperiment::colData(se)$condition == cond)
      sum(is.na(row[idx])) <= thr
    }, logical(1)))
  })

  se[pass, ]
}


#' Filter proteins with too many missing values per condition (fraction-based)
#'
#' Fraction-based variant of filter_missing(). A protein is retained if at
#' least one condition has a proportion of valid (non-NA) replicates >=
#' min_valid_fraction. This threshold scales automatically with replicate count,
#' unlike the absolute `thr` parameter.
#'
#' Example equivalences:
#'   min_valid_fraction = 0.75, n = 4 replicates → requires >= 3 valid
#'   min_valid_fraction = 0.75, n = 5 replicates → requires >= 4 valid
#'
#' @param se                  SummarizedExperiment with colData$condition
#' @param min_valid_fraction  Minimum fraction of valid replicates required in
#'                            the best-covered condition (default 0.75)
#' @return Filtered SummarizedExperiment
filter_missing_frac <- function(se, min_valid_fraction = 0.75) {
  mat <- SummarizedExperiment::assay(se, "intensity")
  conditions <- unique(SummarizedExperiment::colData(se)$condition)

  pass <- apply(mat, 1, function(row) {
    any(vapply(conditions, function(cond) {
      idx <- which(SummarizedExperiment::colData(se)$condition == cond)
      mean(!is.na(row[idx])) >= min_valid_fraction
    }, logical(1)))
  })

  se[pass, ]
}


#' Barplot: proteins detected per sample
#'
#' Replaces DEP::plot_numbers().
#'
#' @param se SummarizedExperiment
#' @return ggplot
plot_sample_numbers <- function(se) {
  mat <- SummarizedExperiment::assay(se, "intensity")
  df <- data.frame(
    sample     = colnames(mat),
    n_detected = colSums(!is.na(mat)),
    condition  = SummarizedExperiment::colData(se)$condition
  )

  ggplot2::ggplot(df, ggplot2::aes(x = sample, y = n_detected, fill = condition)) +
    ggplot2::geom_col(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = n_detected), vjust = -0.3, size = 3) +
    ggplot2::labs(
      title = "Proteins detected per sample",
      x     = "Sample",
      y     = "Number of proteins"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' Cumulative coverage: fraction of proteins detected in at least N samples
#'
#' Replaces DEP::plot_coverage().
#'
#' @param se SummarizedExperiment
#' @return ggplot
plot_protein_coverage <- function(se) {
  mat <- SummarizedExperiment::assay(se, "intensity")
  n_samps <- ncol(mat)
  n_det <- rowSums(!is.na(mat))

  df <- data.frame(
    n_samples     = seq_len(n_samps),
    frac_proteins = vapply(seq_len(n_samps), function(n) mean(n_det >= n), numeric(1))
  )

  ggplot2::ggplot(df, ggplot2::aes(x = n_samples, y = frac_proteins)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(size = 2.5, color = "steelblue") +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::scale_x_continuous(breaks = seq_len(n_samps)) +
    ggplot2::labs(
      title = "Protein detection coverage",
      x     = "Detected in at least N samples",
      y     = "Fraction of proteins"
    ) +
    ggplot2::theme_minimal()
}


#' Heatmap of missing values across proteins and samples
#'
#' Replaces DEP::plot_missval(). Shows only proteins with at least one NA.
#' Requires ComplexHeatmap.
#'
#' @param se   SummarizedExperiment
#' @param name Legend/heatmap name (must be unique when combining with %v% or +; default "Detected")
#' @return ComplexHeatmap Heatmap object (use draw() to render)
plot_missing_heatmap <- function(se, name = "Detected") {
  mat <- SummarizedExperiment::assay(se, "intensity")
  has_na <- rowSums(is.na(mat)) > 0

  # Binary detection matrix: 1 = detected, 0 = missing
  bin_mat <- (!is.na(mat[has_na, ])) * 1L

  # Colour-code the sample annotation by condition
  conditions <- unique(SummarizedExperiment::colData(se)$condition)
  cond_colors <- setNames(
    RColorBrewer::brewer.pal(max(3, length(conditions)), "Set2")[seq_along(conditions)],
    conditions
  )
  col_annot <- ComplexHeatmap::HeatmapAnnotation(
    condition = SummarizedExperiment::colData(se)$condition,
    col       = list(condition = cond_colors)
  )

  ComplexHeatmap::Heatmap(
    bin_mat,
    name = name,
    col = c("0" = "white", "1" = "steelblue"),
    top_annotation = col_annot,
    show_row_names = FALSE,
    cluster_columns = FALSE,
    column_title = sprintf(
      "Missing value pattern — %d proteins\n with ≥1 NA (of %d total)",
      sum(has_na), nrow(mat)
    )
  )
}


#' Intensity density: proteins with vs. without missing values
#'
#' Replaces DEP::plot_detect(). If missing values are MNAR (left-censored),
#' proteins that have NAs should show systematically lower intensities.
#'
#' @param se SummarizedExperiment
#' @return ggplot
plot_intensity_detect <- function(se) {
  mat <- SummarizedExperiment::assay(se, "intensity")
  has_na <- rowSums(is.na(mat)) > 0

  df <- data.frame(
    mean_intensity = c(
      rowMeans(mat[has_na, ], na.rm = TRUE),
      rowMeans(mat[!has_na, ], na.rm = TRUE)
    ),
    group = rep(
      c("With missing values", "Fully detected"),
      c(sum(has_na), sum(!has_na))
    )
  )

  ggplot2::ggplot(df, ggplot2::aes(x = mean_intensity, fill = group, color = group)) +
    ggplot2::geom_density(alpha = 0.35) +
    ggplot2::labs(
      title = "Intensity distribution: detected vs. partially missing",
      x = "Mean log2 intensity",
      y = "Density",
      fill = NULL, color = NULL
    ) +
    ggplot2::theme_minimal()
}


# ── 3. NORMALIZATION — PRONE ──────────────────────────────────────────────────
#
# Install PRONE once before using these functions:
#   renv::install("daisybio/PRONE")
#   renv::snapshot()
#
# PRONE evaluates multiple normalization methods on the same SE and exposes
# diagnostic plots for method selection. Supported methods (as of PRONE 1.x):
#   "GlobalMean", "GlobalMedian", "Quantile", "VSN",
#   "LoessF", "LoessC", "RLR", "EigenMS", "NOMAD", "MBQN"

#' Run all (or selected) PRONE normalization methods on a filtered SE
#'
#' @param se      SummarizedExperiment filtered by filter_missing(); assay
#'                named "intensity" in log2 scale
#' @param methods Character vector of PRONE method names; NULL = all available
#' @return PRONE SummarizedExperiment (one assay per method, plus the raw assay)
normalize_prone <- function(se, methods = NULL) {
  # PRONE expects the assay to be named "log2_intensity" (verify against PRONE docs)
  if (!"log2_intensity" %in% SummarizedExperiment::assayNames(se)) {
    SummarizedExperiment::assayNames(se)[
      SummarizedExperiment::assayNames(se) == "intensity"
    ] <- "log2_intensity"
  }

  if (is.null(methods)) methods <- PRONE::get_normalization_methods()

  # Run normalization — verify the exact function name against PRONE documentation
  # Expected signature: PRONE::normalize(se, methods = methods)
  se_norm <- PRONE::normalize_se(se, methods = methods)
  se_norm
}


#' Diagnostic plots for PRONE normalization comparison
#'
#' Generates density, PCA, intra-group CV, and log-ratio distribution plots
#' for all normalization methods stored in a PRONE SE.
#'
#' @param se_norm PRONE SummarizedExperiment (output of normalize_prone())
#' @return Named list of ggplot objects
plot_normalization_prone <- function(se_norm) {
  # Verify exact function names against PRONE documentation
  list(
    densities  = PRONE::plot_densities(se_norm),
    pca        = PRONE::plot_pca(se_norm),
    cv         = PRONE::plot_intragroup_CV(se_norm),
    log_ratio  = PRONE::plot_log_ratio_distribution(se_norm)
  )
}


#' Extract one normalization method from a PRONE SE as a plain SE
#'
#' @param se_norm PRONE SummarizedExperiment
#' @param method  Name of the normalization method to extract (e.g. "VSN")
#' @return SummarizedExperiment with a single assay named "intensity"
select_normalization <- function(se_norm, method) {
  available <- SummarizedExperiment::assayNames(se_norm)
  if (!method %in% available) {
    stop("Method '", method, "' not found. Available: ", paste(available, collapse = ", "))
  }

  se_out <- se_norm
  # PRONE stores each assay as a data.table, not a base matrix. Coerce back to
  # a matrix (with dimnames restored from se_norm) so downstream code can rely
  # on plain matrix semantics (e.g. mat[, idx, drop = FALSE]).
  norm_mat <- as.matrix(SummarizedExperiment::assay(se_norm, method))
  dimnames(norm_mat) <- dimnames(se_norm)
  # Move the chosen assay into the "intensity" slot and drop the rest
  SummarizedExperiment::assay(se_out, "intensity") <- norm_mat
  SummarizedExperiment::assays(se_out) <-
    SummarizedExperiment::assays(se_out)["intensity"]
  # Add metadata to indicate which normalization was selected
  S4Vectors::metadata(se_out)$selected_normalization <- method
  se_out
}


# ── 4. QC — POST-NORMALIZATION ────────────────────────────────────────────────

#' PCA plot from a SummarizedExperiment
#'
#' Replaces DEP::analyze_dep() + DEP::plot_pca(). Only complete-case proteins
#' (no NAs in any sample) are used for the PCA. Use n_top to restrict to the
#' most variable proteins when the full matrix is large.
#'
#' @param se        SummarizedExperiment
#' @param group_col colData column to colour points by (default "condition")
#' @param n_top     Use top N most variable proteins; NULL = all complete cases
#' @param pc_x      PC index for the x-axis (default 1)
#' @param pc_y      PC index for the y-axis (default 2)
#' @param title     Plot title
#' @return ggplot
plot_pca <- function(se,
                     group_col = "condition",
                     n_top = NULL,
                     pc_x = 1,
                     pc_y = 2,
                     title = "PCA",
                     subtitle = NULL) {
  mat <- SummarizedExperiment::assay(se, "intensity")
  mat <- mat[complete.cases(mat), ]

  # Optionally restrict to the most variable proteins
  if (!is.null(n_top) && n_top < nrow(mat)) {
    row_vars <- apply(mat, 1, var)
    mat <- mat[order(row_vars, decreasing = TRUE)[seq_len(n_top)], ]
  }

  pca <- stats::prcomp(t(mat), scale. = TRUE, center = TRUE)
  var_exp <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  df <- as.data.frame(pca$x[, c(pc_x, pc_y)])
  colnames(df) <- c("PCx", "PCy")
  df$group <- SummarizedExperiment::colData(se)[[group_col]]
  df$sample <- rownames(df)
  df$replicate <- as.factor(SummarizedExperiment::colData(se)$replicate)

  # scale_shape_discrete() (the ggplot2 default) tops out at 6 shapes and
  # silently drops points beyond that. The pch 15-25 range packs several
  # near-duplicate shapes (16/19/20/21 all read as "circle", 15/22 as
  # "square"), so instead we reuse ggplot2's own default 6-shape set (circle,
  # triangle, square, plus, square-cross, asterisk) and extend it with two
  # more visually distinct families (diamond, inverted triangle). Circle,
  # triangle, square and diamond are filled directly by pch; plus/square-cross
  # /asterisk have no interior to fill; the inverted triangle (25) is
  # fillable and rendered solid via fill = after_scale(colour) below.
  n_reps <- nlevels(df$replicate)
  shape_values <- c(16, 17, 15, 3, 7, 8, 18, 25)

  ggplot2::ggplot(df, ggplot2::aes(
    x = PCx, y = PCy,
    color = group, shape = replicate,
    label = sample
  )) +
    ggplot2::geom_point(
      ggplot2::aes(fill = ggplot2::after_scale(colour)),
      size = 4, alpha = 0.85, stroke = 0.3
    ) +
    ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
    ggplot2::scale_shape_manual(values = shape_values[seq_len(n_reps)]) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x        = paste0("PC", pc_x, " (", var_exp[pc_x], "%)"),
      y        = paste0("PC", pc_y, " (", var_exp[pc_y], "%)"),
      color    = group_col
    ) +
    ggplot2::theme_minimal()
}


#' Mean-SD plot: SD vs. rank of mean intensity
#'
#' Replaces DEP::meanSdPlot(). A normalization that removes the mean-variance
#' trend should show a flat lowess line across the intensity range.
#'
#' @param se    SummarizedExperiment
#' @param title Plot title
#' @return ggplot
plot_mean_sd <- function(se, title = "Mean-SD plot") {
  mat <- SummarizedExperiment::assay(se, "intensity")
  means <- rowMeans(mat, na.rm = TRUE)
  sds <- apply(mat, 1, sd, na.rm = TRUE)

  df <- data.frame(rank_mean = rank(means), sd = sds)

  ggplot2::ggplot(df, ggplot2::aes(x = rank_mean, y = sd)) +
    ggplot2::geom_hex(bins = 60) +
    ggplot2::geom_smooth(method = "loess", color = "red", se = FALSE) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::labs(title = title, x = "Rank of mean intensity", y = "SD") +
    ggplot2::theme_minimal()
}


# ── 7. IMPUTATION QC ─────────────────────────────────────────────────────────

#' Intensity density before and after each imputation method
#'
#' Replaces DEP::plot_imputation(). Each panel shows the density of
#' intensities for one method alongside the pre-imputation distribution.
#'
#' @param se_norm      SummarizedExperiment before imputation (has NAs)
#' @param imputed_list Named list of post-imputation SummarizedExperiment objects
#' @return ggplot
plot_imputation_distributions <- function(se_norm, imputed_list) {
  collect_long <- function(se, label) {
    mat <- SummarizedExperiment::assay(se, "intensity")
    data.frame(
      intensity = as.vector(mat),
      sample    = rep(colnames(mat), each = nrow(mat)),
      method    = label
    )
  }

  df <- dplyr::bind_rows(
    collect_long(se_norm, "Before imputation"),
    dplyr::bind_rows(
      lapply(names(imputed_list), function(nm) collect_long(imputed_list[[nm]], nm))
    )
  )
  df <- df[!is.na(df$intensity), ]

  ggplot2::ggplot(df, ggplot2::aes(x = intensity, color = sample)) +
    ggplot2::geom_density(alpha = 0.6, linewidth = 0.5) +
    ggplot2::facet_wrap(~method, ncol = 2) +
    ggplot2::labs(
      title = "Intensity distribution before and after imputation",
      x     = "log2 intensity",
      y     = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
}


#' SD before vs. SD after imputation scatter (per protein per condition)
#'
#' A good imputation preserves per-protein SD (ideal: slope ≈ 1, intercept ≈ 0).
#' Methods that inflate or deflate SD will distort the downstream limma model.
#'
#' @param se_norm      SummarizedExperiment before imputation (has NAs)
#' @param imputed_list Named list of post-imputation SummarizedExperiment objects
#' @return list:
#'   $plot    ggplot scatter faceted by method with lm fit and identity line
#'   $metrics data.frame with slope, intercept, R², mean delta, RMSE per method
plot_sd_before_after <- function(se_norm, imputed_list) {
  sd_before <- .compute_sd_per_protein_condition(se_norm) %>%
    dplyr::rename(sd_before = sd)

  sd_after_all <- dplyr::bind_rows(lapply(names(imputed_list), function(nm) {
    .compute_sd_per_protein_condition(imputed_list[[nm]]) %>%
      dplyr::rename(sd_after = sd) %>%
      dplyr::mutate(method = nm)
  }))

  df <- dplyr::left_join(sd_after_all, sd_before, by = c("name", "condition"))

  # Per-method linear model: sd_after ~ sd_before
  metrics <- df %>%
    dplyr::filter(!is.na(sd_before), !is.na(sd_after)) %>%
    dplyr::group_by(method) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      model  = purrr::map(data, ~ lm(sd_after ~ sd_before, data = .x)),
      glance = purrr::map(model, broom::glance),
      tidy   = purrr::map(model, broom::tidy)
    ) %>%
    tidyr::unnest(tidy) %>%
    dplyr::filter(term %in% c("(Intercept)", "sd_before")) %>%
    tidyr::unnest_wider(glance, names_sep = "_") %>%
    dplyr::select(method, term, estimate, glance_r.squared) %>%
    tidyr::pivot_wider(names_from = term, values_from = estimate) %>%
    dplyr::rename(
      intercept = `(Intercept)`,
      slope     = sd_before,
      r_sq      = glance_r.squared
    ) %>%
    # Add mean delta and RMSE for a fuller picture of distortion
    dplyr::left_join(
      df %>%
        dplyr::mutate(delta = sd_after - sd_before) %>%
        dplyr::group_by(method) %>%
        dplyr::summarize(
          mean_delta = mean(delta, na.rm = TRUE),
          rmse       = sqrt(mean(delta^2, na.rm = TRUE)),
          .groups    = "drop"
        ),
      by = "method"
    ) %>%
    dplyr::mutate(
      label = sprintf(
        "slope=%.2f  int=%.2f\nR²=%.2f  RMSE=%.2f",
        slope, intercept, r_sq, rmse
      )
    )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = sd_before, y = sd_after)) +
    ggplot2::geom_point(alpha = 0.4, size = 0.8) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      color = "red", linetype = "dashed"
    ) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
    ggplot2::geom_text(
      data = metrics,
      ggplot2::aes(x = Inf, y = -Inf, label = label),
      hjust = 1.05, vjust = -0.2,
      size = 2.8, inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~method, scales = "free") +
    ggplot2::labs(
      title = "SD before vs. after imputation  (red = identity, blue = lm fit)",
      x     = "SD before imputation",
      y     = "SD after imputation"
    ) +
    ggplot2::theme_minimal()

  list(plot = p, metrics = as.data.frame(dplyr::select(metrics, -label)))
}


#' Stacked bar: fraction of proteins classified as MNAR / MAR / complete per condition
#'
#' @param mnar_table data.frame from classify_mnar_mar()
#' @return ggplot
plot_mnar_distribution <- function(mnar_table) {
  df <- mnar_table %>%
    dplyr::group_by(condition) %>%
    dplyr::summarize(
      n_total    = dplyr::n(),
      n_mnar     = sum(MNAR_flag),
      n_mar      = sum(!MNAR_flag & num_NAs > 0),
      n_complete = sum(num_NAs == 0),
      .groups    = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols      = c(n_mnar, n_mar, n_complete),
      names_to  = "category",
      values_to = "count"
    ) %>%
    dplyr::mutate(
      category = factor(
        category,
        levels = c("n_mnar", "n_mar", "n_complete"),
        labels = c("MNAR", "MAR", "Complete")
      ),
      fraction = count / n_total
    )

  ggplot2::ggplot(df, ggplot2::aes(x = condition, y = fraction, fill = category)) +
    ggplot2::geom_col(position = "stack", color = "white") +
    ggplot2::geom_text(
      ggplot2::aes(label = count),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::scale_fill_manual(
      values = c(MNAR = "#E74C3C", MAR = "#F39C12", Complete = "#2ECC71")
    ) +
    ggplot2::labs(
      title = "Missing value classification per condition",
      x     = "Condition",
      y     = "Fraction of proteins",
      fill  = NULL
    ) +
    ggplot2::theme_minimal()
}


# Internal helper: per-protein per-condition SD (NA-aware)
.compute_sd_per_protein_condition <- function(se) {
  mat <- SummarizedExperiment::assay(se, "intensity")
  conditions <- unique(SummarizedExperiment::colData(se)$condition)

  dplyr::bind_rows(lapply(conditions, function(cond) {
    idx <- which(SummarizedExperiment::colData(se)$condition == cond)
    sub <- mat[, idx, drop = FALSE]
    data.frame(
      name = rownames(mat),
      condition = cond,
      sd = apply(sub, 1, sd, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}


# ── 8. PCA — OUTLIER DETECTION ────────────────────────────────────────────────

#' Robust PCA outlier detection using PCAHubert (rrcov)
#'
#' Classical PCA is sensitive to outliers; PCAHubert minimises their influence
#' by fitting the subspace to the majority of the data. Samples with high
#' Score Distance (SD) or Orthogonal Distance (OD) relative to the bulk are
#' flagged in a DD-plot. Review before finalising the DE analysis.
#'
#' Requires a fully imputed SE (no NAs in the assay).
#'
#' @param se SummarizedExperiment (fully imputed)
#' @param k  Number of robust PCs to compute (default 5)
#' @return list:
#'   $plot     ggplot DD-plot (score vs. orthogonal distance)
#'   $outliers data.frame with Sample, Condition, Score_Distance,
#'             Orthogonal_Distance, Outlier columns
run_pca_hubert <- function(se, k = 5,
                           title = "Robust PCA outlier detection \nPCAHubert",
                           subtitle = "Samples in red deviate \nfrom the main cluster") {
  mat <- t(SummarizedExperiment::assay(se, "intensity")) # samples × proteins

  fit <- rrcov::PcaHubert(mat, k = k)

  df <- data.frame(
    Sample              = rownames(mat),
    Condition           = SummarizedExperiment::colData(se)$condition,
    Score_Distance      = fit@sd,
    Orthogonal_Distance = fit@od,
    Outlier             = !fit@flag,
    stringsAsFactors    = FALSE
  )

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Score_Distance, y = Orthogonal_Distance,
      color = Outlier, label = Sample
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "steelblue", "TRUE" = "#E74C3C"),
      labels = c("FALSE" = "Normal", "TRUE" = "Outlier")
    ) +
    ggplot2::labs(
      title    = title,
      subtitle = subtitle,
      x        = "Score Distance",
      y        = "Orthogonal Distance",
      color    = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 7),
      plot.subtitle = ggplot2::element_text(size = 6),
      axis.title = ggplot2::element_text(size = 6),
      legend.text = ggplot2::element_text(size = 6),
    )

  list(plot = p, outliers = df)
}


#' PCA grid: one panel per imputation method
#'
#' Calls plot_pca() on each element of imputed_list and combines the results
#' into a patchwork grid. Useful for visually selecting the best method.
#'
#' @param imputed_list Named list of post-imputation SummarizedExperiment objects
#' @param group_col    colData column to colour points by (default "condition")
#' @return patchwork ggplot
plot_pca_comparison <- function(imputed_list, group_col = "condition") {
  labels <- paste0(LETTERS[seq_along(imputed_list)], ". ", names(imputed_list))

  plots <- Map(function(se, lbl) {
    plot_pca(se, group_col = group_col, title = lbl)
  }, imputed_list, labels)

  patchwork::wrap_plots(plots, ncol = 3) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = "PCA comparison across imputation methods"
    ) &
    ggplot2::theme(legend.position = "bottom")
}
