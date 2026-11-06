#!/usr/bin/env Rscript
# Build a small test dataset from the PXD031582 ProteinGroups table (CLN3 interactome).
# Data are NOT versioned: the output goes to testdata/ (git-ignored) and is derived from
# rawdata/ProteinGroups_Project_PXD031582.xlsx (download it from PRIDE, project PXD031582).
#
# Usage: Rscript bin/make_test_data.R [input.xlsx] [output.tsv] [n_proteins]
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[1] else "rawdata/ProteinGroups_Project_PXD031582.xlsx"
output <- if (length(args) >= 2) args[2] else "testdata/proteinGroups_PXD031582_subset.tsv"
n_keep <- if (length(args) >= 3) as.integer(args[3]) else 600L
set.seed(1)

x <- readxl::read_xlsx(input)
lfq <- grep("^LFQ intensity ", names(x), value = TRUE)
stopifnot(length(lfq) == 12)
m <- as.matrix(x[, lfq]); m[is.na(m)] <- 0
n_valid <- rowSums(m > 0)

# Mix of well-quantified proteins and proteins with missing values, so that both
# the MAR and MNAR imputation branches are exercised.
full <- which(n_valid >= 10)
partial <- which(n_valid >= 4 & n_valid < 10)
pick <- c(sample(full, min(length(full), ceiling(n_keep * 0.7))),
          sample(partial, min(length(partial), floor(n_keep * 0.3))))
out <- x[sort(pick), c("accession", "genename", "Description", lfq)]
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(out, output)
message(sprintf("Wrote %d proteins x %d LFQ columns to %s", nrow(out), length(lfq), output))
