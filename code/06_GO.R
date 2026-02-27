# libraries ####

source("../code/00_packages.R")
# source("../code/global_variables.R")


# Load data ####
# Load dep_analysis object
load(paste0(output_path, "RData/dep_analysis.RData"))
# Load data results tables
data_results <- readxl::read_xlsx(path = paste0(output_path, "tables/data_results.xlsx"))


# Gene Ontology ####
# Generate a list with UP and DOWN genes for each comparison
# Use comparisons vector as input


if (independent_UPDOWN == 1) {
  UPDOWN_dataframes <- lapply(comparisons, function(comp) {
    # Extract 'UP' data
    name_UP <- paste0("names_", comp, "_UP")
    df_UP <- data.frame(ID = (data_results %>% dplyr::filter(!!sym(paste0(comp, "_diffexpressed")) == "UP"))$ID)
    # Extract 'DOWN' data
    name_DN <- paste0("names_", comp, "_DN")
    df_DN <- data.frame(ID = (data_results %>% dplyr::filter(!!sym(paste0(comp, "_diffexpressed")) == "DOWN"))$ID)
    # Assign names to data frames
    setNames(list(df_UP, df_DN), c(name_UP, name_DN))
  })
} else {
  UPDOWN_dataframes <- lapply(comparisons, function(comp) {
    # Extract 'UP' and 'DOWN' data
    name_UPDOWN <- paste0("names_", comp, "_UPDOWN")
    df_UPDOWN <- data.frame(ID = (data_results %>%
      dplyr::filter(!!sym(paste0(comp, "_diffexpressed")) %in% c("UP", "DOWN")))$ID)
    # Assign name to data frame
    setNames(list(df_UPDOWN), name_UPDOWN)
  })
}

# Flatten the list of lists into a single list
UPDOWN_dataframes_name <- unlist(UPDOWN_dataframes, recursive = FALSE) # Remove one level in the list
UPDOWN_dataframes <- unlist(UPDOWN_dataframes_name, recursive = FALSE) # Repeat function twice for remove one level more of the list for GO

# GO functions ####

### EnrichGO functions ####

# GO terms list
# Function to enrichGO for the dataframes of a list
perform_enrichGO <- function(gene_set) {
  results <- clusterProfiler::enrichGO(
    gene = gene_set,
    OrgDb = organism,
    keyType = "UNIPROT",
    ont = ont,
    pAdjustMethod = "none",
    readable = TRUE
  )
  return(results)
}

# Apply the enrichGO function to each gene set in differential_genes

assign(paste0("results_list_enrichGO_", ont), lapply(UPDOWN_dataframes, perform_enrichGO))
assign(paste0("results_list_enrichGO_", ont, "_names"), names(paste0("results_list_enrichGO_", ont)))

save(list = paste0("results_list_enrichGO_", ont), file = paste0(output_path, "RData/results_list_enrichGO_", ont, ".RData"))


# load results_list_go_ALL RData file
# load(paste0(output_path,"RData/results_list_enrichGO_",ont,".RData"))

results_list_GO <- get(paste0("results_list_enrichGO_", ont))

# GO Barplots function
perform_barplots <- function(x) {
  if (length(x@result$Count) == 0) { # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    if (ont == "ALL") {
      results <- barplot(x, showCategory = 30) + facet_grid(ONTOLOGY ~ ., scale = "free")
    } else {
      results <- barplot(x, showCategory = 30)
    }
    return(results)
  }
}

# GO Dotplots function
perform_dotplots <- function(x) {
  if (length(x@result$Count) == 0) { # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    if (ont == "ALL") {
      results <- enrichplot::dotplot(x, showCategory = 30) + facet_grid(ONTOLOGY ~ ., scale = "free")
    } else {
      results <- enrichplot::dotplot(x, showCategory = 30)
    }
    return(results)
  }
}

# GO cnetplots function
perform_cnetplots <- function(x) {
  if (length(x@result$Count) == 0) { # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    ## remove redundent GO terms
    x2 <- simplify(x)
    results <- enrichplot::cnetplot(x2)
    # For circular Gene Concept map
    # results <- enrichplot::cnetplot(x2, circular = TRUE, colorEdge = TRUE)
    return(results)
  }
}

# UpSet plot function
perform_upsetplot <- function(x) {
  if (length(x@result$Count) == 0) { # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    results <- enrichplot::upsetplot(x)
    return(results)
  }
}
# image_number <-129

# HeatMap plot function
perform_heatmapplot <- function(x) {
  if (length(x@result$Count) == 0) { # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    results <- enrichplot::heatplot(x)
    return(results)
  }
}


# Tree plot function
perform_treeplot <- function(go_result) {
  if (is.null(go_result) || nrow(go_result@result) == 0) {
    return(NULL)
  }
  edox <- pairwise_termsim(go_result)
  # tree1 <- treeplot(edox)
  tree2 <- treeplot(edox) # nota del warning
  tree2
  # aplot::plot_list(tree1, tree2, tag_levels = "A")
}


# image_number <- 150

### Barplot ####
barplot_results <- lapply(results_list_GO, perform_barplots)
# Iterate over dotplot_results and print each plot individually
for (i in seq_along(barplot_results)) {
  if (is.null(barplot_results[[i]])) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_GO)[i]) # Remove "names_"
    plot_name_hyphen <- gsub("\\.ID", "", plot_name) # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen) # Replace "_" with " "

    p <- barplot_results[[i]] +
      ggplot2::ggtitle(paste0("Bar plot for ", plot_name))
    print(p)
    filename <- paste0(output_path, "figures/enrichGO/", image_number + i, "_barplot_", plot_name_hyphen) # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename, pdf_extension), p, width = 8) # Vectorial format
    ggsave(paste0(filename, tiff_extension), p, width = 8) # Tiff format
  }
}
image_number <- image_number + i

### Dotplot ####

dotplot_results <- lapply(results_list_GO, perform_dotplots)
# Iterate over dotplot_results and print each plot individually
for (i in seq_along(dotplot_results)) {
  if (is.null(dotplot_results[[i]])) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_GO)[i]) # Remove "names_"
    plot_name_hyphen <- gsub("\\.ID", "", plot_name) # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen) # Replace "_" with " "

    p <- dotplot_results[[i]] +
      ggplot2::ggtitle(paste0("Dot plot for ", plot_name))
    print(p)
    filename <- paste0(output_path, "figures/enrichGO/", image_number + i, "_Dotplot_", plot_name_hyphen) # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename, pdf_extension), p, width = 8) # Vectorial format
    ggsave(paste0(filename, tiff_extension), p, width = 8) # Tiff format
  }
}
image_number <- image_number + i

### CNET plot ####

cnetplot_results <- lapply(results_list_GO, perform_cnetplots)

for (i in seq_along(cnetplot_results)) {
  if (is.null(cnetplot_results[[i]])) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_GO)[i]) # Remove "names_"
    plot_name_hyphen <- gsub("\\.ID", "", plot_name) # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen) # Replace "_" with " "

    p <- cnetplot_results[[i]] +
      ggplot2::ggtitle(paste0("CNET plot for ", plot_name))
    print(p)
    filename <- paste0(output_path, "figures/enrichGO/", image_number + i, "_cnet_", plot_name_hyphen) # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename, pdf_extension), p, width = 8, height = 8, units = "in") # Vectorial format
    ggsave(paste0(filename, tiff_extension), p, width = 8, height = 8, units = "in") # Tiff format
  }
}
image_number <- image_number + i


### UPSET plot ####
upsetplot_results <- lapply(results_list_GO, perform_upsetplot)

for (i in seq_along(upsetplot_results)) {
  if (is.null(upsetplot_results[[i]])) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_GO)[i]) # Remove "names_"
    plot_name_hyphen <- gsub("\\.ID", "", plot_name) # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen) # Replace "_" with " "

    p <- upsetplot_results[[i]] +
      ggplot2::ggtitle(paste0("UPSET plot for ", plot_name)) +
      geom_bar(aes(stat = "identity", width = 0.2))
    print(p)
    filename <- paste0(output_path, "figures/enrichGO/", image_number + i, "_upset_", plot_name_hyphen) # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename, pdf_extension), p, width = 8) # Vectorial format
    ggsave(paste0(filename, tiff_extension), p, width = 8) # Tiff format
  }
}
image_number <- image_number + i


### Heatmap plot ####
heatmapplot_results <- lapply(results_list_GO, perform_heatmapplot)

for (i in seq_along(heatmapplot_results)) {
  if (is.null(heatmapplot_results[[i]])) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_GO)[i]) # Remove "names_"
    plot_name_hyphen <- gsub("\\.name$", "", plot_name) # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen) # Replace "_" with " "

    p <- heatmapplot_results[[i]] +
      ggplot2::ggtitle(paste0("HeatMap plot for ", plot_name))
    print(p)
    filename <- paste0(output_path, "figures/enrichGO/", image_number + i, "_heatmap_", plot_name_hyphen) # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename, pdf_extension), p, width = 8, dpi = 300) # Vectorial format
    ggsave(paste0(filename, tiff_extension), p, width = 8, dpi = 300) # Tiff format
  }
}
image_number <- image_number + i


## Tree plots ####
treeplot_results <- lapply(results_list_GO, perform_treeplot)

for (i in seq_along(treeplot_results)) {
  p <- treeplot_results[[i]]
  if (is.null(p)) {
    cat("❌ No se generó plot para: ", names(results_list_GO)[i], "\n")
    next
  }

  # Limpiar el nombre del archivo
  plot_name <- names(results_list_GO)[i] %>%
    gsub("^names_", "", .) %>%
    gsub("\\.name$", "", .)

  filename <- paste0(output_path, "figures/enrichGO/", sprintf("%03d", image_number + i), "_treeplot_", plot_name)

  # Añadir título y exportar
  p <- p + ggtitle(paste("TreePlot for", gsub("_", " ", plot_name)))
  print(p)
  ggsave(paste0(filename, pdf_extension), p, width = 8, height = 6, dpi = 300)
  ggsave(paste0(filename, tiff_extension), p, width = 8, height = 6, dpi = 300)
}
image_number <- image_number + i

# Look at this page follow it.
# https://bioc.ism.ac.jp/packages/3.7/bioc/vignettes/enrichplot/inst/doc/enrichplot.html#bar-plot

## Lolliplots ####
# Make dataframe list from GO results #
# Retrieve the results data frame list


results_df_GO <- map(results_list_GO, ~ .x@result)

results_df_GO %>%
  write_xlsx(paste0(output_path, "tables/results_df_enrichGO_ALL.xlsx"))

results_df_GO <- lapply(results_df_GO, function(df) {
  df <- df %>%
    mutate(GeneRatio = strsplit(GeneRatio, "/") %>%
      map_dbl(~ as.numeric(.x[1]) / as.numeric(.x[2])))
  return(df)
})

lolliplot <- function(data_name, df_list, file_prefix = NULL) {
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

  plot_name <- gsub("^names_", "", names(df_list)[non_empty_indices[data_name]]) # Remove "names_"
  plot_name <- gsub("_vs_", " vs ", plot_name) # Replace "_vs_" with " vs "
  plot_name <- gsub("_", " ", plot_name) # Replace "_" with " "
  plot_name <- gsub("UP$", "UP", plot_name) # Remove "UP" from the end

  plot <- df %>%
    dplyr::mutate(Description = reorder(Description, GeneRatio)) %>%
    top_n(15, GeneRatio) %>%
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

  if (!is.null(file_prefix)) {
    tiff_filename <- paste0(file_prefix, image_number + i, "_Lolliplot", "_", gsub(" ", "_", plot_name), ".tiff")
    pdf_filename <- paste0(file_prefix, image_number + i, "_Lolliplot", "_", gsub(" ", "_", plot_name), ".pdf")
    ggsave(filename = tiff_filename, plot = plot, device = "tiff", width = 8, dpi = 300)
    ggsave(filename = pdf_filename, plot = plot, device = "pdf", width = 8, dpi = 300)
  }
  print(plot)
  return(plot)
}

# Apply lolliplot function
for (i in seq_along(results_df_GO)) {
  lolliplot_result <- lolliplot(i, results_df_GO, file_prefix = paste0(output_path, "figures/enrichGO/"))
}

image_number <- image_number + i
