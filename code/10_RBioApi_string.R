

# libraries ####
source("../code/00_packages.R")
source("../code/global_variables.R")

# Load data ####
#Load dep_analysis object
load(paste0(output_path,"RData/dep_analysis.RData"))
#Load data results tables
data_results <- readxl::read_xlsx(path = paste0(output_path,"tables/data_results.xlsx"))

# Gene String Network ####
# Generate a list with UP and DOWN genes for each comparison
# Use comparisons vector as input
UPDOWN_dataframes <- lapply(comparisons, function(comp) {
  # Extract 'UP' data
  name_UP <- paste0("names_", comp, "_UP")
  df_UP <- data.frame(name = (data_results %>% dplyr::filter(!!sym(paste0(comp, "_diffexpressed")) == "UP"))$name)
  # Extract 'DOWN' data
  name_DN <- paste0("names_", comp, "_DN")
  df_DN <- data.frame(name = (data_results %>% dplyr::filter(!!sym(paste0(comp, "_diffexpressed")) == "DOWN"))$name)
  # Assign names to data frames
  setNames(list(df_UP, df_DN), c(name_UP, name_DN))
})

# Flatten the list of lists into a single list
UPDOWN_dataframes_name <- unlist(UPDOWN_dataframes, recursive = FALSE)# Remove one level in the list
# Export results ####
save(UPDOWN_dataframes_name, file = paste0(output_path,"RData/UPDOWN_dataframes_name.RData")) #File for Strings DB

# Load data ####
#Load dep_analysis object
load(paste0(output_path,"RData/UPDOWN_dataframes_name.RData"))

#species: 9606 for Human, 10090 for mouse, 7955 for zebrafish
perform_string_map <- function(gene_set) {
  proteins_mapped <- rbioapi::rba_string_map_ids(ids = gene_set$name,
                                         species = 7955)
  graph_1 <- rba_string_network_image(ids = proteins_mapped$stringId,
                                      image_format = "highres_image",
                                      species = 7955,
                                      save_image = paste0(output_path,"figures/string/"),
                                      required_score = 500,
                                      network_flavor = "evidence",
                                      network_type = "physical")
  
  return(graph_1)
}
# Apply the enrichGO function to each gene set in differential_genes
results_list_string_ALL <- lapply(UPDOWN_dataframes_name, perform_string_map)


#Save results
save(results_list_string_ALL, file = paste0(output_path,"RData/results_list_string_ALL.RData"))

# load results_list_go_ALL RData file
load(paste0(output_path,"RData/results_list_string_ALL.RData"))


