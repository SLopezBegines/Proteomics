
# Load data ####
#Load dep_analysis object
load(paste0(output_path,"RData/dep_analysis.RData"))
#Load data results tables
data_results <- readxl::read_xlsx(path = paste0(output_path,"tables/data_results.xlsx"))

## Generate list of genes ####

KEGG_input_dataframes <- lapply(comparisons, function(comp) {
  # Extract 'UP' and 'DOWN' data
  name_UPDOWN <- paste0("names_", comp, "_adj")
  df_UPDOWN <- data.frame(name = sig_adjusted_data %>% dplyr::filter(get(paste0(comp, "_diffexpressed")) %in% c("UP", "DOWN")) %>% pull(name),
                          ID = sig_adjusted_data %>% dplyr::filter(get(paste0(comp, "_diffexpressed")) %in% c("UP", "DOWN")) %>% pull(ID),
                          foldchange = sig_adjusted_data %>% dplyr::filter(get(paste0(comp, "_diffexpressed")) %in% c("UP", "DOWN")) %>% pull(!!sym(paste0(comp, "_ratio"))),
                          pvalue = sig_adjusted_data %>% dplyr::filter(get(paste0(comp, "_diffexpressed")) %in% c("UP", "DOWN")) %>% pull(!!sym(paste0(comp, "_p.val"))))
  
  
  # Assign name to data frame
  setNames(list(df_UPDOWN), name_UPDOWN)
})




KEGG_input_dataframes <- unlist(KEGG_input_dataframes, recursive = FALSE)# Keep it for StringDB


# Define a function to process each element in the KEGG_input_dataframes list
process_gene_list <- function(df) {
  # Extract foldchange and ID columns
  original_gene_list <- df$foldchange
  names(original_gene_list) <- df$ID
  
  # Remove NA values
  gene_list <- na.omit(original_gene_list)
  
  # Sort the gene_list in decreasing order
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(gene_list)
}
# Apply the function to each element in the KEGG_input_dataframes list
GSE_gene_lists <- lapply(KEGG_input_dataframes, process_gene_list)

# gse pathway ####
## gseGO function ####
# Function to gseGO for the dataframes of a list
perform_gseGO <- function(gene_set, type) {
  results <- clusterProfiler::gseGO(geneList = gene_set,
                                    ont = type,
                                    keyType = "UNIPROT", 
                                    #nPerm = 10000, 
                                    minGSSize = 3, 
                                    maxGSSize = 800, 
                                    pvalueCutoff = 0.05, 
                                    verbose = TRUE, 
                                    OrgDb = organism, 
                                    pAdjustMethod = "none")
  return(results)
}

# Apply the enrichGO function to each gene set in differential_genes
results_list_gseGO_ALL_adj <- lapply(GSE_gene_lists, perform_gseGO, type= "ALL")
#Save results
save(results_list_gseGO_ALL_adj, file = paste0(output_path,"RData/results_list_gseGO_ALL_adj.RData"))

# load results_list_go_ALL RData file
load(paste0(output_path,"RData/results_list_gseGO_ALL_adj.RData"))

## Plots ####

# GO Dotplots function
perform_dotplots <- function(x) {
  if (nrow(x@result) == 0) {  # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    #results <- enrichplot::dotplot(x, showCategory = 30, split=".sign") + facet_grid(ONTOLOGY ~ ., scale = "free")+ scale_color_viridis()
    results <- enrichplot::dotplot(x, showCategory = 10, split=".sign") + 
      facet_grid(.~.sign) + 
      scale_color_viridis() 
    #results <- enrichplot::dotplot(x, showCategory = 30)
    return(results)
  }
}

# GO cnetplots function
perform_cnetplots <- function(x){
  ## remove redundent GO terms
  x2 <- simplify(x)
  #results <- enrichplot::cnetplot(x2)
  #For circular Gene Concept map
  #results <- enrichplot::cnetplot(x2, circular = TRUE, colorEdge = TRUE)
  results <- enrichplot::cnetplot(x2, foldChange=x@geneList, circular = TRUE, colorEdge = TRUE)
  return(results)
}

#HeatMap plot function
perform_heatmapplot <- function(x){
  if (nrow(x@result) == 0) {  # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    results <- enrichplot::heatplot(x, foldChange=x@geneList)
    return(results)
  }
}
# RidgePlots
perform_ridgeplots <- function(x) {
  if (nrow(x@result) == 0) {  # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    #results <- enrichplot::dotplot(x, showCategory = 30, split=".sign") + facet_grid(ONTOLOGY ~ ., scale = "free")+ scale_color_viridis()
    results <- enrichplot::ridgeplot(x) +  
      scale_fill_viridis()
    #results <- enrichplot::dotplot(x, showCategory = 30)
    return(results)
  }
}
# PMCPlots
perform_pmcplots <- function(x) {
  if (nrow(x@result) == 0) {  # Check if x[result][Count] is empty
    cat("No elements found for dataframe.\n")
    return(NULL)
  } else {
    #results <- enrichplot::dotplot(x, showCategory = 30, split=".sign") + facet_grid(ONTOLOGY ~ ., scale = "free")+ scale_color_viridis()
    terms <- x@result$Description[1:3]
    results <- enrichplot::pmcplot(terms, 2010:2023, proportion=FALSE)
    #results <- enrichplot::dotplot(x, showCategory = 30)
    return(results)
  }
}

#image_number = 103

### Dotplot ####

dotplot_results <- lapply(results_list_gseGO_ALL_adj, perform_dotplots)
# Iterate over dotplot_results and print each plot individually
for (i in seq_along(dotplot_results)) {
  if (length(dotplot_results[[i]]) == 0) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_gseGO_ALL_adj)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.name$", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- dotplot_results[[i]] +
      ggplot2::ggtitle(paste0("Dot plot for ", plot_name))
    print(p)
    filename <- paste0(output_path,"figures/gseGO_adj/",image_number +i,"_GSE_Dotplot_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 8, dpi= 300, units = "in")  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 8, dpi= 300, units = "in")  # Tiff format
  }
}
image_number <- image_number +i

### CNET plot ####
'
cnetplot_results <- lapply(results_list_gseGO_ALL, perform_cnetplots)

for (i in seq_along(cnetplot_results)) {
  if (is.list(cnetplot_results[[i]]) == FALSE) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_gseGO_ALL)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.name$", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- cnetplot_results[[i]]+
      ggplot2::ggtitle(paste0("CNET plot for ", plot_name))
    print(p)
    filename <- paste0(output_path,"figures/gseGO/",image_number+i,"_GSE_cnet_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 8, dpi= 300, units = "in")  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 8, dpi= 300, units = "in")  # Tiff format
  }
}
image_number <- image_number+i

'


### Heatmap plot ####
heatmapplot_results <- lapply(results_list_gseGO_ALL_adj, perform_heatmapplot)

for (i in seq_along(heatmapplot_results)) {
  if (is.list(heatmapplot_results[[i]]) == FALSE) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_gseGO_ALL_adj)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.name$", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- heatmapplot_results[[i]]+
      ggplot2::ggtitle(paste0("HeatMap plot for ", plot_name))
    print(p)
    filename <- paste0(output_path,"figures/gseGO_adj/",image_number+i,"_GSE_heatmap_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 8, units = "in", dpi = 300)  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 8, units = "in", dpi = 300)  # Tiff format
  }
}
image_number <- image_number+i

### Ridge plot ####
ridgeplot_results <- lapply(results_list_gseGO_ALL_adj, perform_ridgeplots)

for (i in seq_along(ridgeplot_results)) {
  if (is.list(ridgeplot_results[[i]]) == FALSE) {
    cat("Plot not generated for element ", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_gseGO_ALL_adj)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.name$", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- ridgeplot_results[[i]]+
      ggplot2::ggtitle(paste0("Ridge plot for ", plot_name))
    print(p)
    filename <- paste0(output_path,"figures/gseGO_adj/",image_number+i,"_GSE_ridgeplot_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 8, units = "in", dpi = 300)  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 8, units = "in", dpi = 300)  # Tiff format
  }
}
image_number <- image_number+i

### PMCplot ####
# Number of articles in pubmed with the input ges descriptions
pmcplot_results <- lapply(results_list_gseGO_ALL_adj, perform_pmcplots)

for (i in seq_along(pmcplot_results)) {
  if (is.list(pmcplot_results[[i]]) == FALSE) {
    cat("Plot not generated for element", i, ".\n")
  } else {
    # Extract the desired part of the name
    plot_name <- gsub("^names_", "", names(results_list_gseGO_ALL_adj)[i])  # Remove "names_"
    plot_name_hyphen <- gsub("\\.name$", "", plot_name)  # Remove ".name"
    plot_name <- gsub("_", " ", plot_name_hyphen)  # Replace "_" with " "
    
    p <- pmcplot_results[[i]]+
      ggplot2::ggtitle(paste0("PMC plot for ", plot_name))
    print(p)
    filename <- paste0(output_path,"figures/gseGO_adj/",image_number+i,"_GSE_pmcplot_", plot_name_hyphen)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
    ggsave(paste0(filename,pdf_extension), p, width = 8, height = 8, units = "in", dpi = 300)  # Vectorial format
    ggsave(paste0(filename,tiff_extension), p, width = 8, height = 8, units = "in", dpi = 300)  # Tiff format
  }
}
image_number <- image_number+i

## Lolliplots ####

# Make dataframe list from GO results #
results_df_gseGO_ALL_adj <- map(results_list_gseGO_ALL_adj, ~ .x@result)

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
  
  plot_name <- gsub("^names_", "", names(df_list)[non_empty_indices[data_name]])  # Remove "names_"
  plot_name <- gsub("_vs_", " vs ", plot_name)  # Replace "_vs_" with " vs "
  plot_name <- gsub("_", " ", plot_name)  # Replace "_" with " "
  plot_name <- gsub("UP$", "UP", plot_name)  # Remove "UP" from the end
  
  plot <- df %>%
    dplyr::mutate(Description = reorder(Description, enrichmentScore)) %>%
    top_n(15, enrichmentScore) %>%
    ggplot(aes(x = Description,
               y = enrichmentScore,
               colour = p.adjust)) +
    geom_segment(aes(x = Description,
                     xend = Description,
                     y = 0,
                     yend = enrichmentScore)) +
    geom_point(aes(size = setSize-1), show.legend = TRUE) +
    scale_color_viridis_c(option = "viridis", direction = 1) +
    facet_wrap(ONTOLOGY ~ ., scale = "free", ncol=1)+
    coord_flip() +
    theme_minimal() +
    labs(size = "N. of genes",
         x = "GO term",
         y = "enrichmentScore") + 
    ggtitle(paste("Plot for", plot_name))
  
  if (!is.null(file_prefix)) {
    tiff_filename <- paste0(file_prefix, image_number+i,"_Lolliplot","_", gsub(" ", "_", plot_name), ".tiff")
    pdf_filename <- paste0(file_prefix, image_number+i,"_Lolliplot","_", gsub(" ", "_", plot_name), ".pdf")
    ggsave(filename = tiff_filename, plot = plot, device = "tiff")
    ggsave(filename = pdf_filename, plot = plot, device = "pdf")
  }
  print(plot)
  return(plot) 
}

# Apply lolliplot function
for (i in seq_along(results_df_gseGO_ALL_adj)) {
  lolliplot_result <- lolliplot(i, results_df_gseGO_ALL_adj, file_prefix = paste0(output_path,"figures/gseGO_adj/"))
}



# BP ####
# Apply the enrichGO function to each gene set in differential_genes
results_list_gseGO_BP <- lapply(GSE_gene_lists, perform_gseGO, type= "BP")
#Save results
save(results_list_gseGO_BP, file = paste0(output_path,"RData/results_list_gseGO_BP_adj.RData"))


# load results_list_go_BP RData file
load(paste0(output_path,"RData/results_list_gseGO_BP_adj.RData"))

# CC ####
# Apply the enrichGO function to each gene set in differential_genes
results_list_gseGO_CC <- lapply(GSE_gene_lists, perform_gseGO, type= "CC")
#Save results
save(results_list_gseGO_CC, file = paste0(output_path,"RData/results_list_gseGO_CC_adj.RData"))

# load results_list_go_CC RData file
load(paste0(output_path,"RData/results_list_gseGO_CC_adj.RData"))
# MF ####
# Apply the enrichGO function to each gene set in differential_genes
results_list_gseGO_MF <- lapply(GSE_gene_lists, perform_gseGO, type= "MF")
#Save results
save(results_list_gseGO_MF, file = paste0(output_path,"RData/results_list_gseGO_MF_adj.RData"))

# load results_list_go_MF RData file
load(paste0(output_path,"RData/results_list_gseGO_MF_adj.RData"))

