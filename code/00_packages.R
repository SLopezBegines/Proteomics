# load libraries ####
# install.packages("renv")
# Function to load libraries from CRAN
install_and_load_library <- function(lib_names) {
  missing_libs <- lib_names[!sapply(lib_names, requireNamespace, quietly = TRUE)]

  if (length(missing_libs) > 0) {
    pak::pkg_install(missing_libs, dependencies = TRUE)
  }

  # Load all libraries
  for (lib in lib_names) {
    library(lib, character.only = TRUE)
  }
}
# Function to load libraries from Bioconductor
install_and_load_library_bioconductor <- function(lib_names) {
  missing_libs <- lib_names[!sapply(lib_names, requireNamespace, quietly = TRUE)]

  if (length(missing_libs) > 0) {
    # Install without dependencies on subsequent calls
    BiocManager::install(missing_libs, dependencies = TRUE, update = FALSE, ask = FALSE, force = FALSE)
  }
  # Load all libraries
  for (lib in lib_names) {
    library(lib, character.only = TRUE)
  }
}
# install.packages("BiocManger")

# Usage example
libraries <- c(
  "tidyverse", "writexl", "BiocManager", "ggrepel", "knitr",
  "kableExtra", "openxlsx", "ggfortify", "ggpubr", "rbioapi",
  "mdatools", "pheatmap", "geneset", "rstatix", "genekitr", "rrcov", "venn", "ggpolypath",
  "ggvenn", "ComplexUpset", "eulerr", "ggVennDiagram", "RColorBrewer", "patchwork"
)

libraries_bioconductor <- c(
  "biomaRt", "org.Dr.eg.db", "org.Mm.eg.db", "org.Hs.eg.db", "clusterProfiler", "enrichplot", "tidySummarizedExperiment",
  "DEP", "AnnotationDbi", "topGO", "STRINGdb", "AnnotationHub", "rrvgo",
  "europepmc", "Rgraphviz", "pathview", "limma", "ComplexHeatmap", "DOSE", "viridis", "edgeR", "vsn"
)

install_and_load_library(libraries)
install_and_load_library_bioconductor(libraries_bioconductor)


rm("install_and_load_library_bioconductor", "install_and_load_library", "libraries", "libraries_bioconductor")
# library(org.Dr.eg.db)#For Zebrafish. Change Dr by Mm or Hs for mouse or human
# BiocManager::install("PANTHER.db")
# library(PANTHER.db)
# BiocManager::install("RDAVIDWebService")
# source("https://bioconductor.org/biocLite.R")
# biocLite("RDAVIDWebService")
# library(RDAVIDWebService)
# library(DescTools)
# library(PerformanceAnalytics)
# library(ggforce)
