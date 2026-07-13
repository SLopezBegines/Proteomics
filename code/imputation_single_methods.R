# imputation_single_methods.R
#
# Drop-in replacements for the four single-method DEP/MSnbase::impute() calls,
# using the ORIGINAL engines directly (no DEP, no MSnbase):
#   MinProb, QRILC  -> imputeLCMD   (CRAN)
#   knn             -> impute        (Bioconductor)
#   man             -> reimplemented (DEP::manual_impute logic, per-column median)
#
# Dependencies:
#   install.packages("imputeLCMD")          # pulls in norm, tmvtnorm, etc.
#   BiocManager::install("impute")
#
# All functions:
#   - take a SummarizedExperiment with a log2 assay named "intensity"
#   - return a SummarizedExperiment with that assay imputed
#   - are pure (no disk writes, no global assignments)
#
# Reproducibility: MinProb, man and QRILC are stochastic. Call set.seed()
# immediately before each, or fix a seed argument, to get stable results.


# ── MinProb ───────────────────────────────────────────────────────────────────
# Reproduces: impute(data_norm, fun = "MinProb", q = 0.01)
# Per-sample random draws from a Gaussian centred on the q-quantile of the
# observed values (left-censored / MNAR). Imputes ALL missing values.
impute_minprob <- function(se, q = 0.01, tune.sigma = 1) {
  mat <- as.matrix(SummarizedExperiment::assay(se, "intensity"))
  imp <- imputeLCMD::impute.MinProb(mat, q = q, tune.sigma = tune.sigma)
  if (is.list(imp)) imp <- imp[[1]] # defensive: some versions wrap
  dimnames(imp) <- dimnames(mat)
  out <- se
  SummarizedExperiment::assay(out, "intensity") <- imp
  out
}


# ── QRILC ─────────────────────────────────────────────────────────────────────
# Reproduces: impute(data_norm, fun = "QRILC")
# Quantile Regression Imputation of Left-Censored data (MNAR).
# impute.QRILC returns a LIST; the imputed matrix is element [[1]].
impute_qrilc <- function(se, tune.sigma = 1) {
  mat <- as.matrix(SummarizedExperiment::assay(se, "intensity"))
  imp <- imputeLCMD::impute.QRILC(mat, tune.sigma = tune.sigma)[[1]]
  dimnames(imp) <- dimnames(mat)
  out <- se
  SummarizedExperiment::assay(out, "intensity") <- imp
  out
}


# ── kNN ───────────────────────────────────────────────────────────────────────
# Reproduces: impute(data_norm, fun = "knn", rowmax = 0.9)
# Uses the SAME engine as MSnbase (impute::impute.knn), so results are identical
# given the same rng.seed. Genes/proteins in rows, samples in columns.
# NOTE: this is NOT the same as qc_proteomics.R::impute_knn_mar (VIM::kNN).
impute_knn_msn <- function(se, k = 10, rowmax = 0.9, colmax = 0.8,
                           maxp = 1500, rng.seed = 362436069) {
  mat <- as.matrix(SummarizedExperiment::assay(se, "intensity"))
  res <- impute::impute.knn(mat,
    k = k, rowmax = rowmax, colmax = colmax,
    maxp = maxp, rng.seed = rng.seed
  )
  imp <- res$data
  dimnames(imp) <- dimnames(mat)
  out <- se
  SummarizedExperiment::assay(out, "intensity") <- imp
  out
}


# ── man (manual left-shifted Gaussian) ────────────────────────────────────────
# Reproduces: impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)
# Faithful reimplementation of DEP::manual_impute: for each sample (column),
# missing values are drawn from N(mean = median_col - shift * sd_col,
# sd = sd_col * scale), using only the OBSERVED values of that column.
impute_man <- function(se, shift = 1.8, scale = 0.3) {
  mat <- as.matrix(SummarizedExperiment::assay(se, "intensity"))
  for (j in seq_len(ncol(mat))) {
    col <- mat[, j]
    na_idx <- is.na(col)
    n_na <- sum(na_idx)
    if (n_na == 0L) next
    med <- stats::median(col, na.rm = TRUE)
    s <- stats::sd(col, na.rm = TRUE)
    mat[na_idx, j] <- stats::rnorm(n_na, mean = med - shift * s, sd = s * scale)
  }
  out <- se
  SummarizedExperiment::assay(out, "intensity") <- mat
  out
}


# ── Mixed spi ────────────────────────────────────────
#
# Drop-in replacement for impute_mixed_split() in qc_proteomics.R.
#
# Method (pure R, no DEP / no MSnbase):
#   - Classification MNAR/MAR is per protein-AND-condition (classify_mnar_mar()).
#     A protein is MNAR in a condition if its fraction of missing replicates in
#     that condition is >= fraction_NA. The SAME protein can be MAR in one
#     condition and MNAR in another (e.g. a knockout: observed in WT, absent in KO).
#
#   - MAR positions are imputed with kNN on the FULL matrix (all samples), using
#     impute::impute.knn — the same engine as the old MSnbase mar = "knn", and
#     consistent with the single-method wrapper impute_knn_msn().
#
#   - Order: MAR first, MNAR excluded. impute.knn is run with the MNAR cells
#     still NA, so the synthetic MNAR draws never enter the neighbour/distance
#     computation. Only the MAR positions are copied back from the kNN result;
#     the values kNN invents in MNAR cells are discarded.
#
#   - MNAR positions are then filled with draws from a left-truncated Gaussian:
#       mean = global_min - factor_SD_impute * global_sd
#       sd   = factor_SD_scale * global_sd
#       b    = global_min   (upper truncation: no draw exceeds the detection floor)
#     Global (whole-matrix) parameters, matching the old 03_cleaning script.
#
# Requires (in addition to qc_proteomics.R deps):
#   BiocManager::install("impute")     # impute::impute.knn
#   install.packages("truncnorm")      # truncnorm::rtruncnorm
#
#' Classify each protein-per-condition as MNAR or MAR
#'
#' A protein is flagged MNAR in a condition when the fraction of missing
#' replicates in that condition is >= fraction_NA. Used by impute_mixed_split()
#' to route each protein to the appropriate imputation method.
#'
#' @param se          SummarizedExperiment (post-normalization, has NAs)
#' @param fraction_NA Fraction of missing replicates to call MNAR (default 0.6)
#' @return data.frame with columns: name, condition, num_NAs, frac_NA, MNAR_flag
classify_mnar_mar <- function(se, fraction_NA = 0.6) {
  mat <- SummarizedExperiment::assay(se, "intensity")
  conditions <- unique(SummarizedExperiment::colData(se)$condition)

  results <- lapply(conditions, function(cond) {
    idx <- which(SummarizedExperiment::colData(se)$condition == cond)
    sub <- mat[, idx, drop = FALSE]
    data.frame(
      name = rownames(mat),
      condition = cond,
      num_NAs = rowSums(is.na(sub)),
      frac_NA = rowMeans(is.na(sub)),
      MNAR_flag = rowMeans(is.na(sub)) >= fraction_NA,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}

#' Mixed imputation: per-condition MNAR/MAR routing (full-matrix kNN)
#'
#' @param se               SummarizedExperiment (post-normalization, assay
#'                         "intensity" in log2 scale, contains NAs)
#' @param fraction_NA      MNAR threshold per condition (default 0.6): a protein
#'                         is MNAR in a condition if >= fraction_NA of its
#'                         replicates there are missing
#' @param factor_SD_impute MNAR centre offset below global_min, as a fraction of
#'                         the global SD (default 0.05)
#' @param factor_SD_scale  MNAR spread, as a fraction of the global SD (default 0.3)
#' @param k                Number of neighbours for impute.knn (default 10)
#' @param rowmax           impute.knn rowmax fallback threshold (default 0.9)
#' @param colmax           impute.knn colmax error threshold (default 0.8)
#' @param maxp             impute.knn block size (default 1500)
#' @param seed             Optional integer for reproducibility. Governs both the
#'                         kNN rng.seed and the truncated-normal draws. NULL = use
#'                         impute.knn's default seed and the current RNG state.
#' @return SummarizedExperiment with a fully imputed assay "intensity" and the
#'         wide MNAR classification table in metadata(se)$mnar_table
impute_mixed_split <- function(se,
                               fraction_NA = 0.6,
                               factor_SD_impute = 0.05,
                               factor_SD_scale = 0.3,
                               k = 10,
                               rowmax = 0.9,
                               colmax = 0.8,
                               maxp = 1500,
                               seed = 1234) {
  mat <- as.matrix(SummarizedExperiment::assay(se, "intensity"))
  cond_vec <- SummarizedExperiment::colData(se)$condition
  conditions <- unique(cond_vec)

  # Global parameters for the MNAR left-censored distribution
  global_min <- min(mat, na.rm = TRUE)
  global_sd <- stats::sd(mat, na.rm = TRUE)

  # 1. Classify each protein-per-condition as MNAR or MAR
  mnar_table <- classify_mnar_mar(se, fraction_NA = fraction_NA)

  # 2. Build full-matrix position masks.
  #    is_mnar: cell is NA AND its protein is MNAR-flagged in that column's condition.
  #    is_mar : any other NA.
  is_mnar <- matrix(FALSE, nrow(mat), ncol(mat), dimnames = dimnames(mat))
  for (cond in conditions) {
    idx <- which(cond_vec == cond)
    tbl_c <- mnar_table[mnar_table$condition == cond, ]
    flags <- stats::setNames(tbl_c$MNAR_flag, tbl_c$name)[rownames(mat)] # align to row order
    sub_na <- is.na(mat[, idx, drop = FALSE])
    is_mnar[, idx] <- flags & sub_na # vector recycled down each column
  }
  is_mar <- is.na(mat) & !is_mnar

  mat_imputed <- mat

  # 3. MAR — kNN on the FULL matrix, with MNAR cells still NA (excluded from
  #    distances/donors). Copy back ONLY the MAR positions.
  if (any(is_mar)) {
    rng <- if (is.null(seed)) 362436069 else as.integer(seed)
    knn_res <- impute::impute.knn(mat,
      k = k, rowmax = rowmax,
      colmax = colmax, maxp = maxp, rng.seed = rng
    )
    knn_full <- knn_res$data
    dimnames(knn_full) <- dimnames(mat)
    mat_imputed[is_mar] <- knn_full[is_mar]
  }

  # 4. MNAR — left-truncated Gaussian below the global detection floor.
  n_mnar <- sum(is_mnar)
  if (n_mnar > 0L) {
    if (!is.null(seed)) set.seed(seed) # after kNN, so truncnorm is reproducible
    value_center <- global_min - factor_SD_impute * global_sd
    sd_scale <- factor_SD_scale * global_sd
    mat_imputed[is_mnar] <- truncnorm::rtruncnorm(
      n = n_mnar, a = -Inf, b = global_min,
      mean = value_center, sd = sd_scale
    )
  }

  # 5. Assemble the imputed SE and attach the classification table (wide) as metadata
  mnar_wide <- tidyr::pivot_wider(
    mnar_table,
    names_from  = condition,
    values_from = c(frac_NA, MNAR_flag, num_NAs),
    names_sep   = "_"
  )

  se_imputed <- se
  SummarizedExperiment::assay(se_imputed, "intensity") <- mat_imputed
  S4Vectors::metadata(se_imputed)$mnar_table <- mnar_wide
  se_imputed
}

# ── Example: ─────────────────────
# set.seed(1)
# data_imp           <- impute_minprob(data_norm, q = 0.01)
# data_imp_man_gauss <- impute_man(data_norm, shift = 1.8, scale = 0.3)
# data_imp_knn       <- impute_knn_msn(data_norm, rowmax = 0.9)
# data_imp_QRILC     <- impute_qrilc(data_norm)
# data_mixed         <- impute_mixed_split(se,
#                                         fraction_NA = 0.6,
#                                         factor_SD_impute = 0.05,
#                                         factor_SD_scale = 0.3,
#                                         k = 10,
#                                         rowmax = 0.9,
#                                         colmax = 0.8,
#                                         maxp = 1500,
#                                         seed = 1234))
# imputed_list <- list(
#   MinProb = data_imp, man = data_imp_man_gauss,
#   knn = data_imp_knn, QRILC = data_imp_QRILC, Mixed = data_mixed
# )
# plot_imputation_distributions(data_norm, imputed_list)   # from qc_proteomics.R
# plot_sd_before_after(data_norm, imputed_list)            # from qc_proteomics.R
