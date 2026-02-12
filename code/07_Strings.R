



# libraries ####

source("../code/00_packages.R")
#source("../code/global_variables.R")

# Load data ####
#Load dep_analysis object
load(paste0(output_path,"RData/dep_analysis.RData"))
#Load data results tables
data_results <- readxl::read_xlsx(path = paste0(output_path,"tables/data_results.xlsx"))


# Gene Ontology ####
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

UPDOWN_dataframes_name %>% 
  write_xlsx(paste0(output_path,"tables/","STRING_dataframes_name.xlsx"))
# Load data ####
#Load dep_analysis object
load(paste0(output_path,"RData/UPDOWN_dataframes_name.RData"))



# Strings ####
# load string database
#9606 for Human, 10090 for mouse, 7955 for zebrafish
string_db <- STRINGdb$new(version = "11.5", species = species, score_threshold = 200, input_directory="")
class(string_db)


string_function <- function(x) {
  string_example <- string_db$map(x, "ID", removeUnmappedRows = TRUE)
  dimension <- dim(string_example)[1]
  hits <- string_example$STRING_id[1:dimension]
  link <- as.character(string_db$get_link(hits))
  return(link)
}

STRING_plotnames <- lapply(names(UPDOWN_dataframes_name), function(name) {
  link <- string_function(UPDOWN_dataframes_name[[name]])
  data.frame(Name = name, Link = link, stringsAsFactors = FALSE)
})

### Save String plot links ####
STRING_plotnames_df <- do.call(rbind, STRING_plotnames)
STRING_plotnames_table <- STRING_plotnames_df %>%
  kable("html") %>%
  kable_styling()

print(STRING_plotnames_table)

STRING_plotnames_df %>% 
  write_xlsx(paste0(output_path,"tables/STRING_plotnames_table.xlsx"))
