

# libraries ####

source("./code/00_packages.R")
source("code/global_variables.R")

library(enrichR)


# Load data ####
#Load dep_analysis object
load(paste0(output_path,"RData/dep_analysis.RData"))
#Load data results tables
data_results <- readxl::read_xlsx(path = paste0(output_path,"tables/data_results.xlsx"))


# Gene Ontology ####
# Generate a list with UP and DOWN genes for each comparison
# Use comparisons vector as input
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

# Flatten the list of lists into a single list
UPDOWN_dataframes <- unlist(UPDOWN_dataframes, recursive = FALSE)# Remove one level in the list
UPDOWN_dataframes <- unlist(UPDOWN_dataframes, recursive = FALSE)#Repeat function twice for remove one level more of the list for GO

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
#List of all databases in dbs object
if (websiteLive) dbs <- listEnrichrDbs()
head(dbs)
#Select the desired database

database <- c("KEGG_2016","GO_Biological_Process_2023")

# EnrichRGO terms list
# Function to enrichGO for the dataframes of a list
perform_enrichR <- function(gene_set) {
  results <- enrichR::enrichr(gene_set, 
                              databases = database)
  return(results)
}


# Apply the enrichGO function to each gene set in differential_genes


results_list_enrichGO_ALL <- lapply(UPDOWN_dataframes_name, perform_enrichR)
#results_list_enrichGO_ALL_names <- names(results_list_enrichGO_ALL)

save(results_list_enrichR, file = paste0(output_path,"RData/results_list_enrichR.RData"))

# load results_list_go_ALL RData file
load(paste0(output_path,"RData/results_list_enrichR.RData"))

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}
if (websiteLive) dbs <- listEnrichrDbs()
