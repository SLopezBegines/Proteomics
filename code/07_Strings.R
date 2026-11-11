# STRING protein-protein interaction links ####
# Builds the UP / DOWN protein sets of each comparison from `data_results` and
# queries STRING for one network link per set.

source("../code/00_packages.R")

#' STRING network links for the DE protein sets of each comparison
#'
#' @param data_results       Table from run_de_analysis(); read from `output_path` if NULL.
#' @param comparisons        Character vector of comparison names.
#' @param independent_UPDOWN TRUE: separate UP and DOWN sets; FALSE: one combined set.
#' @param diffexpressed_type "qval", "pval" or "padj_holm_pval".
#' @param top_n_FC           NULL = all DE proteins; n = the n with the largest |FC| per set.
#' @param species            NCBI taxonomy ID (9606 human, 10090 mouse, 7955 zebrafish).
#' @param output_path        Project output folder (RData/ and tables/ are written there).
#' @param string_version,score_threshold STRINGdb version and minimum interaction score.
#' @return Data frame (Name, Link), also written to tables/STRING_plotnames_table.xlsx.
run_STRING <- function(data_results = NULL,
                       comparisons,
                       independent_UPDOWN = TRUE,
                       diffexpressed_type = c("qval", "pval", "padj_holm_pval"),
                       top_n_FC = NULL,
                       species = params$species,
                       output_path = params$output_path,
                       string_version = "12.0",
                       score_threshold = 200) {
  diffexpressed_type <- match.arg(diffexpressed_type)
  if (is.null(data_results)) {
    data_results <- readxl::read_xlsx(file.path(output_path, "tables", "data_results.xlsx"))
  }

  # One data frame (column `ID`, UniProt accession) per comparison and direction.
  # Protein groups are stored as "ACC1;ACC2": split them, STRING needs one ID per row.
  UPDOWN_dataframes_name <- prepare_DE_gene_sets(
    data_results = data_results,
    comparisons = comparisons,
    direction = c("UP", "DOWN"),
    independent_UPDOWN = independent_UPDOWN,
    id_cols = "ID",
    flatten = FALSE,
    diffexpressed_type = diffexpressed_type,
    top_n = top_n_FC
  )
  UPDOWN_dataframes_name <- lapply(UPDOWN_dataframes_name, function(df) {
    data.frame(ID = unique(unlist(strsplit(as.character(df$ID), ";", fixed = TRUE))))
  })

  # Export sets (file for STRING DB)
  save(UPDOWN_dataframes_name, file = file.path(output_path, "RData", "UPDOWN_dataframes_name.RData"))
  writexl::write_xlsx(UPDOWN_dataframes_name, file.path(output_path, "tables", "STRING_dataframes_name.xlsx"))

  string_db <- STRINGdb$new(version = string_version, species = species, score_threshold = score_threshold, input_directory = "")

  # STRING network link of a set of IDs (NA if no ID maps to STRING)
  string_link <- function(x) {
    mapped <- string_db$map(x, "ID", removeUnmappedRows = TRUE)
    if (nrow(mapped) == 0) {
      return(NA_character_)
    }
    as.character(string_db$get_link(mapped$STRING_id))
  }

  links <- do.call(rbind, lapply(names(UPDOWN_dataframes_name), function(name) {
    data.frame(Name = name, Link = string_link(UPDOWN_dataframes_name[[name]]), stringsAsFactors = TRUE)
  }))

  print(links %>% kable("html") %>% kable_styling())
  writexl::write_xlsx(links, file.path(output_path, "tables", "STRING_plotnames_table.xlsx"))
  invisible(links)
}
