
# libraries ####

source("../code/00_packages.R")
#source("../code/global_variables.R")


# Load data ####
#Load dep_analysis object
load(paste0(output_path,"RData/dep_analysis.RData"))
#Load data results tables
# Leer datos
data_results <- readxl::read_xlsx(paste0(output_path, "tables/data_results.xlsx"))
#data_results <- readxl::read_xlsx(paste0("mains/",output_path, "tables/data_results.xlsx"))

# Crear lista de genes UP/DOWN

if(independent_UPDOWN == 1){
  
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
}else{
  
  UPDOWN_dataframes <- lapply(comparisons, function(comp) {
    # Extract 'UP' and 'DOWN' data
    name_UPDOWN <- paste0("names_", comp, "_UPDOWN")
    df_UPDOWN <- data.frame(name = (data_results %>% 
                                    dplyr::filter(!!sym(paste0(comp, "_diffexpressed")) %in% c("UP", "DOWN")))$name)
    # Assign name to data frame
    setNames(list(df_UPDOWN), name_UPDOWN)
  })
}


# Aplanar la lista
UPDOWN_dataframes <- unlist(UPDOWN_dataframes, recursive = FALSE)
UPDOWN_dataframes <- unlist(UPDOWN_dataframes, recursive = FALSE)#Repeat function twice for remove one level more of the list for GO

names(UPDOWN_dataframes) <- gsub("^names_(.*?)(_UP|_DN|_UPDOWN)?(\\.name)?$", "\\1\\2", names(UPDOWN_dataframes))
names(UPDOWN_dataframes) <- gsub("\\.name$", "", names(UPDOWN_dataframes))  # limpieza opcional


# Venn Diagrams ####

remove_na_genes <- function(x) x[!is.na(x)]


venn_plot <- function(gene_list, list_name) {
  gene_list <- lapply(gene_list, remove_na_genes)
  
  n_sets <- length(gene_list)
  palette <- brewer.pal(min(n_sets, 8), "Set2")  # máximo 8 colores
  
  p <- ggvenn::ggvenn(
    gene_list,
    fill_color = palette,
    stroke_size = 0.5,
    text_size = 4
  )

  ggsave(paste0(output_path,"VennDiagram/", list_name, "_venn.pdf"), p, width = 7, height = 7)
  ggsave(paste0(output_path,"VennDiagram/", list_name, "_venn.tiff"), p, width = 7, height = 7)
  
  print(p)
  return(p)
}


# Un solo Venn general con todas las comparaciones
venn_plot(UPDOWN_dataframes, "All_Groups")

# Upset plots ####

library(ComplexUpset)

upset_plot <- function(gene_list, list_name) {
  genes <- unique(unlist(gene_list))
  gene_df <- data.frame(gene = genes)
  
  for (set in names(gene_list)) {
    gene_df[[set]] <- gene_df$gene %in% gene_list[[set]]
  }
  
  p <- upset(gene_df, intersect = names(gene_list), name = "Gene Sets", min_size = 1)
  
  ggsave(paste0(output_path,"VennDiagram/", list_name, "_upset.pdf"), p, width = 10, height = 6)
  ggsave(paste0(output_path,"VennDiagram/", list_name, "_upset.tiff"), p, width = 10, height = 6)
  
  print(p)
  return(p)
}

upset_plot(UPDOWN_dataframes, "All_Groups")

# Intersections ####

calculate_all_intersections <- function(gene_lists, list_name) {
  gene_lists <- lapply(gene_lists, remove_na_genes)
  gene_lists <- gene_lists[lengths(gene_lists) > 0]
  
  if (length(gene_lists) == 0) return(data.frame())
  
  intersection_data <- list()
  
  for (n in 1:length(gene_lists)) {
    combinations <- combn(seq_along(gene_lists), n, simplify = FALSE)
    intersections_n <- lapply(combinations, function(comb) {
      genes <- Reduce(intersect, gene_lists[comb])
      list(
        name = paste(names(gene_lists)[comb], collapse = "_AND_"),
        genes = paste(genes, collapse = ", "),
        num_genes = length(genes)
      )
    })
    intersection_data <- c(intersection_data, intersections_n)
  }
  
  data.frame(
    List_Name = list_name,
    Intersection = sapply(intersection_data, `[[`, "name"),
    Genes = sapply(intersection_data, `[[`, "genes"),
    Num_Genes = sapply(intersection_data, `[[`, "num_genes"),
    stringsAsFactors = FALSE
  )
}

final_intersection_df <- calculate_all_intersections(UPDOWN_dataframes, "All_Groups")
print(final_intersection_df)
write.xlsx(final_intersection_df, paste0(output_path,"VennDiagram/intersections_genes.xlsx"), rowNames = FALSE)


