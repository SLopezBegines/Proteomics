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
# Usage example
cran_packages <- c(
  "tidyverse", "writexl", "ggrepel", "knitr",
  "kableExtra", "openxlsx", "ggfortify", "rbioapi",
  "mdatools", "pheatmap", "rstatix", "rrcov",
  "ggvenn", "ComplexUpset", "RColorBrewer", "patchwork",
  "qvalue", "truncnorm", "broom", "gmp", "mixOmics", "msigdbr", "randomForest",
  "enrichR", "readxl", "rmdformats", "VIM"
)
load_pkgs(cran_packages)
# =============================================================================
# Bioconductor packages
# =============================================================================
message("[PACKAGES] Loading Bioconductor packages...")


packages_bioconductor <- c(
  "BiocVersion",
  "biomaRt", "org.Dr.eg.db", "org.Mm.eg.db", "org.Hs.eg.db", "clusterProfiler", "enrichplot", "STRINGdb",
  "pathview", "limma", "ComplexHeatmap", "viridis", "edgeR", "vsn", "DEP"
)

load_pkgs(packages_bioconductor)

# =============================================================================
# GitHub packages
# =============================================================================
message("[PACKAGES] Loading GitHub packages...")
gh_packages <- c(
  "PRONE" = "daisybio/PRONE"
)

for (pkg_name in names(gh_packages)) {
  pkg_spec <- gh_packages[[pkg_name]]
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    message("[INSTALL] Installing from GitHub: ", pkg_spec)
    pak::pkg_install(pkg_spec, ask = FALSE, upgrade = FALSE)
  }
  if (requireNamespace(pkg_name, quietly = TRUE)) {
    suppressPackageStartupMessages(suppressWarnings(
      library(pkg_name, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE)
    ))
  }
}

#   renv::snapshot()

rm(check_sysreqs, load_pkgs, cran_packages, packages_bioconductor, gh_packages)
