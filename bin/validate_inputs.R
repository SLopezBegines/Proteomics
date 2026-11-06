#!/usr/bin/env Rscript
# Validate the protein table, the design and the contrasts; fail early with clear messages.
# Writes design.validated.csv, contrasts.txt and versions.yml.
this_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
source(file.path(this_dir, "pipeline_common.R"))
suppressPackageStartupMessages({ library(readr); library(readxl) })

a <- parse_args(list(contrasts = ""))
for (k in c("protein_groups", "design")) if (is.null(a[[k]])) die("Missing --", k)
if (!file.exists(a$protein_groups)) die("Protein table not found: ", a$protein_groups)
if (!file.exists(a$design)) die("Design file not found: ", a$design)

design <- read_design(a$design)
need <- c("label", "condition")
if (!all(need %in% names(design))) die("Design must have columns: label, condition (found: ", paste(names(design), collapse = ", "), ")")
if (anyNA(design$label) || anyNA(design$condition)) die("Design has empty label/condition values.")
if (anyDuplicated(design$label)) die("Duplicated sample labels in the design: ", paste(design$label[duplicated(design$label)], collapse = ", "))
if (any(grepl("_vs_", design$condition, fixed = TRUE))) die("Condition names must not contain '_vs_' (it is the contrast separator).")
design$condition <- as.character(design$condition)
n_rep <- table(design$condition)
if (any(n_rep < 2)) die("Conditions with fewer than 2 replicates: ", paste(names(n_rep)[n_rep < 2], collapse = ", "))
if (anyDuplicated(paste(design$condition, design$replicate))) die("condition + replicate must be unique per sample.")

# Every design label must be a column of the protein table
header <- names(read_protein_table(a$protein_groups))
missing_labels <- setdiff(design$label, header)
if (length(missing_labels)) {
  die("Design labels not found as columns in ", basename(a$protein_groups), ": ",
      paste(utils::head(missing_labels, 5), collapse = "; "),
      if (length(missing_labels) > 5) " ..." else "",
      "\n  Available LFQ-like columns: ", paste(utils::head(grep("LFQ|Intensity|intensity", header, value = TRUE), 5), collapse = "; "))
}

conditions <- unique(design$condition)
contrasts <- if (nzchar(a$contrasts)) trimws(strsplit(a$contrasts, ",")[[1]]) else
  apply(utils::combn(conditions, 2), 2, function(x) paste(x, collapse = "_vs_"))
for (cn in contrasts) {
  parts <- strsplit(cn, "_vs_", fixed = TRUE)[[1]]
  if (length(parts) != 2 || !all(parts %in% conditions)) {
    die("Invalid contrast '", cn, "'. Use <condition1>_vs_<condition2> with conditions from the design: ", paste(conditions, collapse = ", "))
  }
}

readr::write_csv(design, "design.validated.csv")
writeLines(contrasts, "contrasts.txt")
cat(sprintf("Inputs OK: %d samples, %d conditions (%s), %d contrasts (%s)\n",
            nrow(design), length(conditions), paste(conditions, collapse = ", "), length(contrasts), paste(contrasts, collapse = ", ")))
writeLines(c('"VALIDATE_INPUTS":', paste0("    r-base: ", paste(R.version$major, R.version$minor, sep = ".")),
             paste0("    readxl: ", as.character(packageVersion("readxl"))), paste0("    readr: ", as.character(packageVersion("readr")))), "versions.yml")
