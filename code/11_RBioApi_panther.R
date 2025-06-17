



# libraries ####
source("../code/00_packages.R")
#source("../code/global_variables.R")

# Load data ####
#Load dep_analysis object
load(paste0(output_path,"RData/dep_analysis.RData"))
#Load data results tables
data_results <- readxl::read_xlsx(path = paste0(output_path,"tables/data_results.xlsx"))


# Generate a list with UP and DOWN genes for each comparison
# Use comparisons vector as input


if(independent_UPDOWN == 1){
  
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
}else{
  
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
UPDOWN_dataframes_name <- unlist(UPDOWN_dataframes, recursive = FALSE)# Remove one level in the list
# Export results ####
save(UPDOWN_dataframes_name, file = paste0(output_path,"RData/UPDOWN_dataframes_name.RData")) #File for Strings DB

# Load data ####
#Load dep_analysis object
load(paste0(output_path,"RData/UPDOWN_dataframes_name.RData"))

# Get list of annotation dataset. Choose one for annot_dataset option into rba_panther_enrich
annots <- rba_panther_info(what = "datasets")


# Enrichment ####
#Statistical overrepresentation test

#species: 9606 for Human, 10090 for mouse, 7955 for zebrafish
perform_enrich <- function(gene_set) {
  enriched <- rba_panther_enrich(genes = gene_set$ID,
                                 organism = species,
                                 annot_dataset = "GO:0008150", #from annots result
                                 cutoff = 0.05)
  
  return(enriched)
}
# Apply the enrichGO function to each gene set in differential_genes
results_panther <- lapply(UPDOWN_dataframes_name, perform_enrich)


#Save results
save(results_panther, file = paste0(output_path,"RData/results_panther.RData"))

# load results_list_go_ALL RData file
load(paste0(output_path,"RData/results_panther.RData"))

# Make dataframe list from Panther results #
GO_panther <- map(results_panther, ~ .x$result)

GO_panther %>% 
  write_xlsx(paste0(output_path,"tables/results_GO_panther.xlsx"))

# Plots ####

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
    dplyr::mutate(term.label = reorder(term.label, fold_enrichment)) %>%
    top_n(15, fold_enrichment) %>%
    ggplot(aes(x = term.label,
               y = fold_enrichment,
               colour = -log10(fdr))) +
    geom_segment(aes(x = term.label,
                     xend = term.label,
                     y = 0,
                     yend = fold_enrichment)) +
    geom_point(aes(size = number_in_list), show.legend = TRUE) +
    scale_color_viridis_c(option = "viridis", direction = 1) +
    coord_flip() +
    theme_minimal() +
    labs(size = "N. of genes",
         x = "GO term",
         y = "Fold enrichment") + 
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
for (i in seq_along(GO_panther)) {
  lolliplot_result <- lolliplot(i, GO_panther, file_prefix = paste0(output_path,"figures/panther/Lolliplot_0",i))
}

'
# Mapping ####

perform_map_panther <- function(gene_set){
  result <- rba_panther_mapping(genes = gene_set$ID,
                                   organism = species)
  return(result)
}

map_panther <- lapply(UPDOWN_dataframes_name, perform_map_panther)


#Save results
save(map_panther, file = "output/RData/map_panther.RData")

# load results_list_go_ALL RData file
load("output/RData/map_panther.RData")

mapped_genes_panther <- map(map_panther, ~ .x$mapped_genes$gene)



# Initialize empty lists to store data
gene_ids <- character()
biological_functions <- character()

# Iterate over the nested list using lapply
lapply(df, function(sublist) {
  # Extract gene ID
  gene <- sublist$accession
  
  # Extract annotation data
  annotations <- sublist$annotation_type_list$annotation_data_type
  
  # Iterate over each annotation data type
  lapply(annotations, function(annotation) {
    # Extract UniProt ID
    uniprot_id <- sub(".*UniProtKB=([^|]+).*", "\\1", gene)
    
    # Extract biological functions
    functions <- annotation$annotation_list$annotation$name
    
    # Append biological functions and UniProt ID to lists
    biological_functions <<- c(biological_functions, functions)
    gene_ids <<- c(gene_ids, rep(uniprot_id, length(functions)))
  })
})

# Create dataframe
result_df <- data.frame(Biological_Function = biological_functions, Gene_ID = gene_ids)

# Rearrange dataframe to have biological functions in the first column
# and all genes mapping to each biological function in the second column
result_df <- aggregate(Gene_ID ~ Biological_Function, data = result_df, paste, collapse = ", ")

# Print dataframe
print(result_df)

################
your_list <- mapped_genes_panther

results_list_go <- list()



your_list <- mapped_genes_panther[[12]]
# Initialize empty lists to store data
gene_ids <- character()
biological_functions <- character()

# Iterate over the nested list using lapply
lapply(your_list, function(sublist) {
  # Extract gene ID
  gene <- sublist$accession
  
  # Extract annotation data
  annotations <- sublist$annotation_type_list$annotation_data_type
  
  # Iterate over each annotation data type
  lapply(annotations, function(annotation) {
    # Extract UniProt ID
    uniprot_id <- sub(".*UniProtKB=([^|]+).*", "\\1", gene)
    
    # Extract biological functions
    functions <- annotation$annotation_list$annotation$name
    
    # Append biological functions and UniProt ID to lists
    biological_functions <<- c(biological_functions, functions)
    gene_ids <<- c(gene_ids, rep(uniprot_id, length(functions)))
  })
})

# Create dataframe
result_df <- data.frame(Biological_Function = biological_functions, Gene_ID = gene_ids)

# Rearrange dataframe to have biological functions in the first column
# and all genes mapping to each biological function in the second column
result_df <- aggregate(Gene_ID ~ Biological_Function, data = result_df, paste, collapse = ", ")

# Print dataframe
#print(result_df)

results_list_go[[12]] <- result_df  # Almacenar los resultados en la lista



##############
# Crear una lista de nombres de las listas contenidas en lista_de_listas
names_list <- names(mapped_genes_panther)

# Iterar sobre los nombres de las listas y asignarlos a los nombres de los dataframes
for (i in 1:length(names_list)) {
  names(results_list_go[i]) <- names_list[i]
}

'

