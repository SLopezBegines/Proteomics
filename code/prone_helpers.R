# PRONE helpers ####
#
# Glue code between the pipeline SummarizedExperiment (assay "intensity",
# linear scale) and PRONE (https://github.com/lisiarend/PRONE).
#
# Contents
#   1. as_prone_se()            adapt a pipeline SE to what PRONE expects
#   2. prone_select_methods()   choose normalization methods valid for a design
#   3. prone_valid_assays()     assays that survive the PCA preconditions
#   4. run_prone_normalization() as_prone_se() + PRONE::normalize_se()


# ── 1. SE adaptation ──────────────────────────────────────────────────────────

#' Adapt a pipeline SummarizedExperiment to PRONE's conventions
#'
#' PRONE needs (i) assays named "raw" (linear scale) and "log2", (ii) assay
#' dimnames-free rows, (iii) a rowData column "Protein.IDs", (iv) a colData
#' column "Column" holding the sample names, and (v) the names of the
#' condition (and optionally batch) colData columns in `metadata(se)`.
#'
#' Row names of the assays are stripped on purpose: PRONE's own outputs are
#' data.tables (no row names) and some helpers (e.g. the zero-variance filter
#' in get_complete_pca_dt()) compare `which(...)` against `integer(0)` with
#' identical(). With row names `which()` returns a *named* integer(0), the guard
#' misfires and the assay is filtered down to 0 rows, which breaks plot_PCA()
#' for that assay only.
#'
#' @param se            SummarizedExperiment with a linear-scale intensity
#'                      assay (zeros already converted to NA).
#' @param assay_name    Character(1). Name of the linear-scale assay in `se`.
#'                      Default "intensity".
#' @param id_col        Character(1). rowData column with the protein
#'                      identifier; copied to rowData "Protein.IDs" (read by
#'                      EigenMS and NormicsVSN/NormicsMedian). Default "ID".
#' @param condition_col Character(1). colData column with the grouping factor
#'                      (read by EigenMS and the intragroup plots). Default
#'                      "condition".
#' @param batch_col     Character(1) or NULL. colData column with the batch.
#'                      NULL (default) = no batch; required only by limBE/IRS.
#' @param refs_col      Character(1) or NULL. colData column flagging
#'                      reference channels; required only by IRS. Default NULL.
#' @return SummarizedExperiment with assays `raw` and `log2`, ready for
#'         PRONE::normalize_se() and the PRONE plotting functions.
as_prone_se <- function(se,
                        assay_name = "intensity",
                        id_col = "ID",
                        condition_col = "condition",
                        batch_col = NULL,
                        refs_col = NULL) {
  stopifnot(
    is.character(assay_name), length(assay_name) == 1L,
    assay_name %in% SummarizedExperiment::assayNames(se),
    id_col %in% colnames(SummarizedExperiment::rowData(se)),
    condition_col %in% colnames(SummarizedExperiment::colData(se))
  )
  for (cc in c(batch_col, refs_col)) {
    if (!cc %in% colnames(SummarizedExperiment::colData(se))) {
      stop("colData column '", cc, "' not found.")
    }
  }

  raw_mat <- as.matrix(SummarizedExperiment::assay(se, assay_name))
  log2_mat <- log2(raw_mat)
  log2_mat[is.infinite(log2_mat)] <- NA
  rownames(raw_mat) <- NULL
  rownames(log2_mat) <- NULL

  se_prone <- se
  SummarizedExperiment::assays(se_prone, withDimnames = FALSE) <- S4Vectors::SimpleList(
    raw  = raw_mat,
    log2 = log2_mat
  )

  SummarizedExperiment::rowData(se_prone)$Protein.IDs <- SummarizedExperiment::rowData(se)[[id_col]]
  SummarizedExperiment::colData(se_prone)$Column <- colnames(se_prone)

  S4Vectors::metadata(se_prone)$condition <- condition_col
  if (!is.null(batch_col)) S4Vectors::metadata(se_prone)$batch <- batch_col
  if (!is.null(refs_col)) S4Vectors::metadata(se_prone)$refs <- refs_col

  se_prone
}


# ── 2. Method selection ───────────────────────────────────────────────────────

#' Normalization methods that are applicable to a given design
#'
#' limBE needs a batch column with >= 2 levels; IRS additionally needs TMT
#' reference channels, so it is never valid for label-free data. Without this
#' filter PRONE::normalize_se() aborts with "No batch specified!".
#'
#' @param se        SummarizedExperiment (pipeline or PRONE one).
#' @param batch_col Character(1) or NULL. colData batch column; NULL or a
#'                  column with a single level drops limBE.
#' @param has_refs  Logical(1). TRUE only for multi-batch TMT designs with
#'                  reference channels (keeps IRS). Default FALSE.
#' @param methods   Character vector of candidate methods; NULL (default) =
#'                  PRONE::get_normalization_methods().
#' @param exclude   Character vector of methods to drop unconditionally.
#' @param verbose   Logical(1). Message which methods were dropped and why.
#' @return Character vector of usable method names.
prone_select_methods <- function(se,
                                 batch_col = NULL,
                                 has_refs = FALSE,
                                 methods = NULL,
                                 exclude = NULL,
                                 verbose = TRUE) {
  if (is.null(methods)) methods <- PRONE::get_normalization_methods()

  n_batches <- if (is.null(batch_col)) {
    0L
  } else {
    dplyr::n_distinct(SummarizedExperiment::colData(se)[[batch_col]])
  }

  drop <- exclude
  if (n_batches < 2L) drop <- c(drop, "limBE")
  if (!has_refs || n_batches < 2L) drop <- c(drop, "IRS")
  drop <- intersect(unique(drop), methods)

  if (verbose && length(drop) > 0) {
    message(
      "[PRONE] Methods dropped: ", paste(drop, collapse = ", "),
      " (batch levels = ", n_batches, ", IRS refs = ", has_refs, ")"
    )
  }
  setdiff(methods, drop)
}


# ── 3. Assays valid for PCA ───────────────────────────────────────────────────

#' Normalized assays that can be passed to PRONE::plot_PCA()
#'
#' Replicates PRONE's PCA preconditions per assay: drop proteins with any NA,
#' drop zero-variance proteins, and require prcomp() to succeed on the rest.
#' Assays that fail (e.g. EigenMS can leave no complete protein) are reported
#' and excluded instead of breaking the whole plot call.
#'
#' @param se_norm  PRONE SummarizedExperiment returned by normalize_se().
#' @param exclude  Character vector of assays never to include. Default
#'                 c("raw", "log2") (the un-normalized inputs).
#' @param verbose  Logical(1). Message the per-assay diagnostics.
#' @return Character vector of assay names where prcomp() succeeded.
prone_valid_assays <- function(se_norm, exclude = c("raw", "log2"), verbose = TRUE) {
  candidates <- setdiff(SummarizedExperiment::assayNames(se_norm), exclude)

  ok <- vapply(candidates, function(a) {
    m <- as.matrix(SummarizedExperiment::assay(se_norm, a))
    m <- m[stats::complete.cases(m), , drop = FALSE]
    zero_var <- apply(m, 1, stats::var) == 0
    m <- m[!zero_var, , drop = FALSE]

    res <- tryCatch(
      {
        stats::prcomp(t(m), scale. = TRUE)
        TRUE
      },
      error = function(e) {
        if (verbose) message("[PRONE] ", a, ": prcomp failed (", conditionMessage(e), ")")
        FALSE
      }
    )
    if (verbose && res) {
      message("[PRONE] ", a, ": ", nrow(m), " complete proteins, ", sum(zero_var), " zero-variance removed")
    }
    res
  }, logical(1))

  candidates[ok]
}


# ── 4. One-call normalization ─────────────────────────────────────────────────

#' Adapt a pipeline SE to PRONE and run the requested normalizations
#'
#' @param se            SummarizedExperiment (see as_prone_se()).
#' @param methods       Character vector of PRONE method names. Pass a single
#'                      method for production runs, or the output of
#'                      prone_select_methods() for an exploratory comparison.
#' @param assay_name,id_col,condition_col,batch_col,refs_col
#'                      Passed to as_prone_se().
#' @return PRONE SummarizedExperiment: assays `raw`, `log2` and one per method.
run_prone_normalization <- function(se,
                                    methods,
                                    assay_name = "intensity",
                                    id_col = "ID",
                                    condition_col = "condition",
                                    batch_col = NULL,
                                    refs_col = NULL) {
  stopifnot(is.character(methods), length(methods) >= 1L)
  unknown <- setdiff(methods, PRONE::get_normalization_methods())
  if (length(unknown) > 0) {
    stop("Unknown PRONE method(s): ", paste(unknown, collapse = ", "))
  }

  se_prone <- as_prone_se(
    se,
    assay_name = assay_name, id_col = id_col, condition_col = condition_col,
    batch_col = batch_col, refs_col = refs_col
  )
  PRONE::normalize_se(se_prone, methods = methods)
}
