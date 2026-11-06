# Shared helpers for the Nextflow stage scripts in bin/.
# These wrap the legacy functions in code/ WITHOUT modifying them: they recreate the global
# environment that mains/*.Rmd sets up (organism, output_path, thresholds, image_number).

# Directory the process was started in (Nextflow task dir). setup_pipeline() changes the
# working directory, so anything that must land here has to use this path.
ORIG_WD <- getwd()

parse_args <- function(defaults = list()) {
  a <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  i <- 1
  while (i <= length(a)) {
    if (startsWith(a[i], "--")) {
      key <- sub("^--", "", a[i])
      if (i < length(a) && !startsWith(a[i + 1], "--")) { out[[key]] <- a[i + 1]; i <- i + 2 }
      else { out[[key]] <- "TRUE"; i <- i + 1 }
    } else i <- i + 1
  }
  out
}

die <- function(...) {
  cat("ERROR:", ..., "\n", file = stderr())
  quit(save = "no", status = 1)
}

organism_table <- list(
  human     = list(kegg = "hsa", taxid = 9606, orgdb = "org.Hs.eg.db"),
  mouse     = list(kegg = "mmu", taxid = 10090, orgdb = "org.Mm.eg.db"),
  zebrafish = list(kegg = "dre", taxid = 7955, orgdb = "org.Dr.eg.db")
)

read_design <- function(path) {
  d <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (!"batch" %in% names(d)) d$batch <- 1
  if (!"replicate" %in% names(d)) d$replicate <- stats::ave(seq_len(nrow(d)), d$condition, FUN = seq_along)
  as.data.frame(d)
}

read_contrasts <- function(path) {
  x <- trimws(readLines(path, warn = FALSE))
  x[nzchar(x)]
}

# Legacy Exp_design layout expected by code/03 and code/04.
build_exp_design <- function(design, dataset) {
  conditions <- unique(as.character(design$condition))
  data.frame(
    label = as.character(design$label),
    condition = factor(design$condition, levels = conditions),
    replicate = as.character(design$replicate),
    batch = design$batch,
    columns_to_rename = paste(design$condition, design$replicate, sep = "_"),
    experiment = dataset,
    stringsAsFactors = FALSE
  )
}

# Read the protein table (xlsx / tsv / csv). Accepts the pre-processed layout used by the
# legacy scripts (accession, genename, Description) and the raw MaxQuant proteinGroups.txt.
read_protein_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  x <- switch(ext,
    xlsx = , xls = readxl::read_excel(path),
    csv = readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    tsv = , txt = readr::read_tsv(path, show_col_types = FALSE, progress = FALSE),
    die("Unsupported input extension '", ext, "' (use .xlsx, .csv, .tsv or .txt)")
  )
  x <- as.data.frame(x, check.names = FALSE)
  is_raw_maxquant <- all(c("Majority protein IDs", "Gene names") %in% names(x)) && !"accession" %in% names(x)
  if (is_raw_maxquant) {
    for (flag in c("Reverse", "Potential contaminant", "Only identified by site")) {
      if (flag %in% names(x)) x <- x[is.na(x[[flag]]) | x[[flag]] != "+", , drop = FALSE]
    }
    x$accession <- sub(";.*$", "", x[["Majority protein IDs"]])
    x$genename <- sub(";.*$", "", x[["Gene names"]])
    if ("Protein names" %in% names(x)) x$Description <- x[["Protein names"]]
  }
  x
}

# Recreate the repository layout the legacy scripts rely on (source("../code/...")), set the
# global variables the Rmd would define, and source the shared helpers in the same order.
setup_pipeline <- function(ctx, state_dir) {
  code_dir <- Sys.getenv("PROTEOMICS_CODE_DIR", "/opt/proteomics/code")
  if (!dir.exists(code_dir)) die("Code directory not found: ", code_dir)
  state_dir <- normalizePath(state_dir, mustWork = TRUE)

  layout <- file.path(tempdir(), "legacy_layout")
  dir.create(file.path(layout, "mains"), recursive = TRUE, showWarnings = FALSE)
  file.symlink(code_dir, file.path(layout, "code"))
  setwd(file.path(layout, "mains"))

  org <- organism_table[[ctx$organism]]
  assign("kegg_organism", org$kegg, envir = .GlobalEnv)
  assign("species", org$taxid, envir = .GlobalEnv)
  assign("organism", org$orgdb, envir = .GlobalEnv)
  assign("output_path", paste0(state_dir, "/"), envir = .GlobalEnv)
  assign("image_number", ctx$image_number, envir = .GlobalEnv)

  source("../code/global_variables.R", echo = FALSE)
  source("../code/aux_functions.R", echo = FALSE)
  source("../code/00_packages.R", echo = FALSE)

  # Pipeline parameters override the defaults of global_variables.R
  assign("p_val", ctx$p_val, envir = .GlobalEnv)
  assign("p_val_low", ctx$p_val_low, envir = .GlobalEnv)
  assign("FC", ctx$fc, envir = .GlobalEnv)
  conds <- levels(ctx$exp_design$condition)
  assign("my_colors", stats::setNames(
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(conds)), conds
  ), envir = .GlobalEnv)
  invisible(NULL)
}

context_file <- function(state_dir) file.path(state_dir, "RData", "pipeline_context.rds")
save_context <- function(ctx, state_dir) {
  ctx$image_number <- get("image_number", envir = .GlobalEnv)
  dir.create(file.path(state_dir, "RData"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ctx, context_file(state_dir))
}
load_context <- function(state_dir) {
  f <- context_file(state_dir)
  if (!file.exists(f)) die("Missing ", f, ": the upstream step did not complete.")
  readRDS(f)
}

write_versions <- function(task, path = file.path(ORIG_WD, "versions.yml")) {
  pv <- function(p) as.character(utils::packageVersion(p))
  lines <- c(
    paste0('"', task, '":'),
    paste0("    r-base: ", paste(R.version$major, R.version$minor, sep = ".")),
    paste0("    bioconductor: ", as.character(BiocManager::version())),
    paste0("    DEP: ", pv("DEP")),
    paste0("    limma: ", pv("limma")),
    paste0("    SummarizedExperiment: ", pv("SummarizedExperiment")),
    paste0("    clusterProfiler: ", pv("clusterProfiler"))
  )
  writeLines(lines, path)
}
