# load libraries ####
# install.packages("renv")
# renv::restore()

# Bootstrap pak itself --------------------------------------------------------
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak", repos = sprintf(
    "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
    .Platform$pkgType, R.Version()$os, R.Version()$arch
  ))
}
library(pak)

# Print system requirements for any missing packages (informational) ----------
check_sysreqs <- function(pkgs) {
  tryCatch(
    {
      reqs <- pak::pkg_sysreqs(pkgs)
      if (length(reqs$packages) > 0) {
        message("[SYSREQS] Missing system packages detected:")
        message("  Run: sudo apt-get install -y ", paste(reqs$packages, collapse = " "))
      }
    },
    error = function(e) invisible(NULL)
  )
}

# Helper: install + load, skipping already-installed -------------------------
# pkgs can use pak prefixes (bioc::Pkg, user/repo) for install;
# library() needs bare names, so strip everything up to and including "::".
load_pkgs <- function(pkgs) {
  bare <- sub("^[^:]+::", "", pkgs)
  missing_idx <- !vapply(bare, requireNamespace, logical(1), quietly = TRUE)
  if (any(missing_idx)) {
    message("[INSTALL] Installing: ", paste(pkgs[missing_idx], collapse = ", "))
    pak::pkg_install(pkgs[missing_idx], ask = FALSE, upgrade = FALSE)
  }
  invisible(lapply(bare, function(p) {
    suppressPackageStartupMessages(suppressWarnings(
      library(p, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
    ))
  }))
}

# =============================================================================
# CRAN packages# =============================================================================
message("[PACKAGES] Loading CRAN packages...")
cran_packages <- c(
  "tidyverse", "writexl", "ggrepel", "knitr",
  "kableExtra", "openxlsx", "ggfortify", "rbioapi",
  "mdatools", "pheatmap", "rstatix", "rrcov",
  "ggvenn", "ComplexUpset", "RColorBrewer", "patchwork",
  "qvalue", "truncnorm", "broom", "gmp", "mixOmics", "msigdbr", "randomForest",
  "enrichR", "readxl", "rmdformats", "VIM",
  "viridis", "imputeLCMD", "circlize"
)
load_pkgs(cran_packages)
# =============================================================================
# Bioconductor packages
# =============================================================================
message("[PACKAGES] Loading Bioconductor packages...")


# bioc:: makes pak resolve against the Bioconductor release matching the
# running R (R 4.6 -> Bioc 3.23). BiocVersion is pulled in as a dependency.
# DEP is intentionally absent: it was removed from Bioc 3.23.
packages_bioconductor <- paste0("bioc::", c(
  "biomaRt", "org.Dr.eg.db", "org.Mm.eg.db", "org.Hs.eg.db",
  "clusterProfiler", "enrichplot", "STRINGdb", "pathview", "limma",
  "ComplexHeatmap", "edgeR", "vsn", "SummarizedExperiment", "impute", "S4Vectors", "PRONE"
))

load_pkgs(packages_bioconductor)

#install.packages("devtools")
#library(devtools)
#devtools::install_github("cit-bioinfo/mMCP-counter")

#   renv::snapshot()

rm(check_sysreqs, load_pkgs, cran_packages, packages_bioconductor)
