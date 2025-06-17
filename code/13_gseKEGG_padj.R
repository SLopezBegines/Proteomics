

# Load data ####
#Load dep_analysis object
load(paste0(output_path,"RData/dep_analysis.RData"))
#Load data results tables
data_results <- readxl::read_xlsx(path = paste0(output_path,"tables/data_results.xlsx"))

#Select only significant protein according p-value.
sig_adjusted_data <- readxl::read_xlsx(path = paste0(output_path,"tables/sig_adjusted_data.xlsx"))


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


# gseKEGG ####
## gseKEGG function ####
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
###########

# Define a function to process each element in the KEGG_input_dataframes list
KEGG_gene_list <- function(df) {
  #Take UNIPROT gene IDs
  ids_ENSEMBL <- df$ID |> 
    bitr(fromType = "UNIPROT", toType = "ENSEMBL", OrgDb=organism)
  # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
  dedup_ids <- ids_ENSEMBL[!duplicated(ids_ENSEMBL$ENSEMBL), ]
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  df2 <- left_join(df, dedup_ids, by = c("ID" = "UNIPROT"))
  df2 <- df2[!is.na(df2$ENSEMBL), ]
  
  # Extract foldchange and ID columns
  original_gene_list <- df2$foldchange
  names(original_gene_list) <- df2$ENSEMBL
  # Remove NA values
  gene_list <- na.omit(original_gene_list)
  # Sort the gene_list in decreasing order
  gene_list <- sort(gene_list, decreasing = TRUE)
  return(gene_list)
}
# Apply the function to each element in the KEGG_input_dataframes list
KEGG_gene_lists <- lapply(KEGG_input_dataframes, KEGG_gene_list)


# Function to gseKEGG for the dataframes of a list
perform_gseKEGG <- function(gene_set) {
  results <- clusterProfiler::gseKEGG(geneList = gene_set,
                                      organism     = kegg_organism,
                                      nPerm        = 10000,
                                      minGSSize    = 3,
                                      maxGSSize    = 800,
                                      pvalueCutoff = 0.05,
                                      pAdjustMethod = "none",
                                      keyType       = "uniprot")
  return(results)
}

# Apply the enrichGO function to each gene set in differential_genes
results_list_gseKEGG_ALL_padj <- lapply(GSE_gene_lists, perform_gseKEGG)
#Save results
save(results_list_gseKEGG_ALL_padj, file = paste0(output_path,"RData/results_list_gseKEGG_ALL_padj.RData"))

# load results_list_go_ALL RData file
load(paste0(output_path,"RData/results_list_gseKEGG_ALL_padj.RData"))


## Lolliplots ####

# Make dataframe list from GO results #
results_df_gseKEGG_ALL <- map(results_list_gseKEGG_ALL, ~ .x@result)
# Write each dataframe in the list to a separate sheet in an Excel file
write_xlsx(results_df_gseKEGG_ALL, path = paste0(output_path,"tables/results_df_gseKEGG_ALL_padj.xlsx"))

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
    #facet_wrap(ONTOLOGY ~ ., scale = "free", ncol=1)+
    coord_flip() +
    theme_minimal() +
    labs(size = "N. of genes",
         x = "GO term",
         y = "enrichmentScore") + 
    ggtitle(paste("Plot for", plot_name))
  
  if (!is.null(file_prefix)) {
    tiff_filename <- paste0(file_prefix, "_", gsub(" ", "_", plot_name), ".tiff")
    pdf_filename <- paste0(file_prefix, "_", gsub(" ", "_", plot_name), ".pdf")
    ggsave(filename = tiff_filename, plot = plot, device = "tiff")
    ggsave(filename = pdf_filename, plot = plot, device = "pdf")
  }
  print(plot)
  return(plot) 
}

# Apply lolliplot function
for (i in seq_along(results_df_gseKEGG_ALL)) {
  lolliplot_result <- lolliplot(i, results_df_gseKEGG_ALL, file_prefix = paste0(output_path,"figures/KEGG_GO_adj/Lolliplot_0",i))
}




#Select manually the desired pathway
'
#KEGG PATHWAY
library(pathview)

# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=GSE_gene_lists$names_CLN3_Lux1_vs_WT_UP, pathway.id="dre01200", species = kegg_organism)

print(dme)
# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=GSE_gene_lists$names_CLN3_Lux1_vs_WT_UP, pathway.id="dre00620", species = kegg_organism, kegg.native = T)
knitr::include_graphics("dre00620.pathview.png")
'