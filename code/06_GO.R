# libraries ####

source("../code/00_packages.R")
# source("../code/global_variables.R")


# GO enrichment function ####
# Runs clusterProfiler::enrichGO() on UP/DOWN/UPDOWN gene sets for each
# comparison. Uses `data_results` from the global environment if available,
# otherwise loads it from `output_path`.

run_GO_enrichment <- function(data_results = NULL,
                              comparisons,
                              direction = c("UP", "DOWN"),
                              independent_UPDOWN = TRUE,
                              diffexpressed_type = c("qval", "pval", "padj_holm_pval"),
                              ont = "ALL",
                              organism,
                              output_path,
                              keyType = "UNIPROT",
                              pAdjustMethod = "none") {
  # Use data_results from the global environment if available, otherwise load it
  if (is.null(data_results)) {
    if (exists("data_results", envir = .GlobalEnv)) {
      data_results <- get("data_results", envir = .GlobalEnv)
    } else {
      data_results <- readxl::read_xlsx(path = paste0(output_path, "tables/data_results.xlsx"))
    }
  }

  # Build UP / DOWN gene sets for each comparison
  gene_sets <- prepare_DE_gene_sets(
    data_results = data_results,
    comparisons = comparisons,
    direction = direction,
    independent_UPDOWN = independent_UPDOWN,
    id_cols = "ID",
    flatten = TRUE,
    diffexpressed_type = diffexpressed_type
  )

  # Run enrichGO for each gene set
  perform_enrichGO <- function(gene_set) {
    clusterProfiler::enrichGO(
      gene = gene_set,
      OrgDb = organism,
      keyType = keyType,
      ont = ont,
      pAdjustMethod = pAdjustMethod,
      readable = TRUE
    )
  }

  results_list_GO <- lapply(gene_sets, perform_enrichGO)

  save(results_list_GO, file = paste0(output_path, "RData/results_list_enrichGO_", ont, ".RData"))

  return(results_list_GO)
}


# GO plots function ####
# Generates bar/dot/cnet/upset/heatmap/tree/lolliplots for a results list
# produced by run_GO_enrichment(), and exports the GO results tables.

plot_GO_results <- function(results_list_GO,
                            direction = c("UP", "DOWN"),
                            ont = "ALL",
                            output_path,
                            output_dir = "enrichGO",
                            showCategory = 30,
                            top_n_terms = 15,
                            width = 8,
                            height = 12) {
  direction <- match.arg(direction, choices = c("UP", "DOWN"), several.ok = TRUE)

  # Keep only the gene sets matching the requested direction(s)
  patterns <- c()
  if ("UP" %in% direction) patterns <- c(patterns, "_UP(\\.|$)", "_UPDOWN(\\.|$)")
  if ("DOWN" %in% direction) patterns <- c(patterns, "_DN(\\.|$)", "_UPDOWN(\\.|$)")
  results_list_GO <- results_list_GO[grepl(paste(unique(patterns), collapse = "|"), names(results_list_GO))]

  fig_dir <- file.path(output_path, "figures", output_dir)

  ### Plot functions ####

  # GO Barplots function
  perform_barplots <- function(x) {
    if (length(x@result$Count) == 0) {
      cat("No elements found for dataframe.\n")
      return(NULL)
    }
    if (ont == "ALL") {
      barplot(x, showCategory = showCategory) + facet_grid(ONTOLOGY ~ ., scale = "free")
    } else {
      barplot(x, showCategory = showCategory)
    }
  }

  # GO Dotplots function
  perform_dotplots <- function(x) {
    if (length(x@result$Count) == 0) {
      cat("No elements found for dataframe.\n")
      return(NULL)
    }
    if (ont == "ALL") {
      enrichplot::dotplot(x, showCategory = showCategory) + facet_grid(ONTOLOGY ~ ., scale = "free")
    } else {
      enrichplot::dotplot(x, showCategory = showCategory)
    }
  }

  # GO cnetplots function
  perform_cnetplots <- function(x) {
    if (length(x@result$Count) == 0) {
      cat("No elements found for dataframe.\n")
      return(NULL)
    }
    ## remove redundant GO terms
    x2 <- simplify(x)
    enrichplot::cnetplot(x2)
  }

  # UpSet plot function
  perform_upsetplot <- function(x) {
    if (length(x@result$Count) == 0) {
      cat("No elements found for dataframe.\n")
      return(NULL)
    }
    enrichplot::upsetplot(x)
  }

  # HeatMap plot function
  perform_heatmapplot <- function(x) {
    if (length(x@result$Count) == 0) {
      cat("No elements found for dataframe.\n")
      return(NULL)
    }
    enrichplot::heatplot(x)
  }

  # Tree plot function
  perform_treeplot <- function(go_result) {
    if (is.null(go_result) || nrow(go_result@result) == 0) {
      return(NULL)
    }
    edox <- pairwise_termsim(go_result)
    treeplot(edox)
  }

  ### Barplot ####
  barplot_results <- lapply(results_list_GO, perform_barplots)
  for (i in seq_along(barplot_results)) {
    if (is.null(barplot_results[[i]])) {
      cat("Plot not generated for element ", i, ".\n")
    } else {
      plot_name_hyphen <- gsub("\\.ID$", "", gsub("^names_", "", names(results_list_GO)[i]))
      plot_name <- gsub("_", " ", plot_name_hyphen)

      p <- barplot_results[[i]] +
        ggplot2::ggtitle(paste0("Bar plot for ", plot_name))
      save_plot(paste0("barplot_", plot_name_hyphen),
        p,
        output_dir = fig_dir,
        width = width,
        height = height
      )
    }
  }

  ### Dotplot ####
  dotplot_results <- lapply(results_list_GO, perform_dotplots)
  for (i in seq_along(dotplot_results)) {
    if (is.null(dotplot_results[[i]])) {
      cat("Plot not generated for element ", i, ".\n")
    } else {
      plot_name_hyphen <- gsub("\\.ID$", "", gsub("^names_", "", names(results_list_GO)[i]))
      plot_name <- gsub("_", " ", plot_name_hyphen)

      p <- dotplot_results[[i]] +
        ggplot2::ggtitle(paste0("Dot plot for ", plot_name))
      save_plot(paste0("dotplot_", plot_name_hyphen),
        p,
        output_dir = fig_dir,
        width = width,
        height = height
      )
    }
  }

  ### CNET plot ####
  cnetplot_results <- lapply(results_list_GO, perform_cnetplots)
  for (i in seq_along(cnetplot_results)) {
    if (is.null(cnetplot_results[[i]])) {
      cat("Plot not generated for element ", i, ".\n")
    } else {
      plot_name_hyphen <- gsub("\\.ID$", "", gsub("^names_", "", names(results_list_GO)[i]))
      plot_name <- gsub("_", " ", plot_name_hyphen)

      p <- cnetplot_results[[i]] +
        ggplot2::ggtitle(paste0("CNET plot for ", plot_name))
      save_plot(paste0("cnet_plot_", plot_name_hyphen),
        p,
        output_dir = fig_dir,
        width = width,
        height = height
      )
    }
  }

  ### UPSET plot ####
  upsetplot_results <- lapply(results_list_GO, perform_upsetplot)
  for (i in seq_along(upsetplot_results)) {
    if (is.null(upsetplot_results[[i]])) {
      cat("Plot not generated for element ", i, ".\n")
    } else {
      plot_name_hyphen <- gsub("\\.ID$", "", gsub("^names_", "", names(results_list_GO)[i]))
      plot_name <- gsub("_", " ", plot_name_hyphen)

      p <- upsetplot_results[[i]] +
        ggplot2::ggtitle(paste0("UPSET plot for ", plot_name)) +
        geom_bar(aes(stat = "identity", width = 0.2))
      save_plot(paste0("upset_", plot_name_hyphen),
        p,
        output_dir = fig_dir,
        width = width,
        height = height
      )
    }
  }

  ### Heatmap plot ####
  heatmapplot_results <- lapply(results_list_GO, perform_heatmapplot)
  for (i in seq_along(heatmapplot_results)) {
    if (is.null(heatmapplot_results[[i]])) {
      cat("Plot not generated for element ", i, ".\n")
    } else {
      plot_name_hyphen <- gsub("\\.ID$", "", gsub("^names_", "", names(results_list_GO)[i]))
      plot_name <- gsub("_", " ", plot_name_hyphen)

      p <- heatmapplot_results[[i]] +
        ggplot2::ggtitle(paste0("HeatMap plot for ", plot_name))
      save_plot(paste0("heatmap_", plot_name_hyphen),
        p,
        output_dir = fig_dir,
        width = width,
        height = height
      )
    }
  }

  ## Tree plots ####
  treeplot_results <- lapply(results_list_GO, perform_treeplot)
  for (i in seq_along(treeplot_results)) {
    p <- treeplot_results[[i]]
    if (is.null(p)) {
      cat("No se generó plot para: ", names(results_list_GO)[i], "\n")
    } else {
      plot_name_hyphen <- gsub("\\.ID$", "", gsub("^names_", "", names(results_list_GO)[i]))
      plot_name <- gsub("_", " ", plot_name_hyphen)

      p <- p + ggtitle(paste("TreePlot for", plot_name))
      save_plot(paste0("treeplot_", plot_name_hyphen),
        p,
        output_dir = fig_dir,
        width = width,
        height = height
      )
    }
  }

  # Look at this page follow it.
  # https://bioc.ism.ac.jp/packages/3.7/bioc/vignettes/enrichplot/inst/doc/enrichplot.html#bar-plot

  ## Lolliplots ####

  # Retrieve the results data frame list and export it
  results_df_GO <- map(results_list_GO, ~ .x@result)

  results_df_GO %>%
    write_xlsx(paste0(output_path, "tables/results_df_enrichGO_", ont, ".xlsx"))

  results_df_GO_ratio <- lapply(results_df_GO, function(df) {
    df %>%
      mutate(GeneRatio = strsplit(GeneRatio, "/") %>%
        map_dbl(~ as.numeric(.x[1]) / as.numeric(.x[2])))
  })

  lolliplot <- function(data_name, df_list) {
    non_empty_indices <- which(sapply(df_list, function(x) nrow(x) > 0))

    if (length(non_empty_indices) == 0) {
      message("All data frames are empty. No plots generated.")
      return(NULL)
    }

    if (length(non_empty_indices) < data_name) {
      message("Selected index is out of range. No plot generated.")
      return(NULL)
    }

    df <- df_list[[non_empty_indices[data_name]]]

    if (nrow(df) == 0) {
      message("Selected data frame is empty. No plot generated.")
      return(NULL)
    }

    plot_name_hyphen <- gsub("\\.ID$", "", gsub("^names_", "", names(df_list)[non_empty_indices[data_name]]))
    plot_name <- gsub("_", " ", plot_name_hyphen)

    plot <- df %>%
      dplyr::mutate(Description = reorder(Description, GeneRatio)) %>%
      top_n(top_n_terms, GeneRatio) %>%
      ggplot(aes(
        x = Description,
        y = GeneRatio,
        colour = p.adjust
      )) +
      geom_segment(aes(
        x = Description,
        xend = Description,
        y = 0,
        yend = GeneRatio
      )) +
      geom_point(aes(size = Count), show.legend = TRUE) +
      scale_color_viridis_c(option = "viridis", direction = 1) +
      facet_wrap(ONTOLOGY ~ ., scale = "free", ncol = 1) +
      coord_flip() +
      theme_minimal() +
      labs(
        size = "N. of genes",
        x = "GO term",
        y = "Gene Ratio"
      ) +
      ggtitle(paste("Plot for", plot_name))

    save_plot(paste0("lolliplot_", plot_name_hyphen),
      plot,
      output_dir = fig_dir,
      width = width,
      height = height
    )

    return(plot)
  }

  # Apply lolliplot function
  for (i in seq_along(results_df_GO_ratio)) {
    lolliplot(i, results_df_GO_ratio)
  }

  invisible(results_df_GO)
}


# Run pipeline ####

# results_list_GO <- run_GO_enrichment(
#  comparisons = comparisons,
#  direction = c("UP", "DOWN"),
#  independent_UPDOWN = (independent_UPDOWN == 1),
#  ont = ont,
#  organism = organism,
#  output_path = output_path
# )

# results_df_GO <- plot_GO_results(
#  results_list_GO = results_list_GO,
#  direction = c("UP", "DOWN"),
#  ont = ont,
#  output_path = output_path
# )
