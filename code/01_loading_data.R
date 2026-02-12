

# libraries ####
source("../code/00_packages.R")

# load data ####


loading_data <- function(prot_data){
  name_df <- deparse(substitute(prot_data))
  prot_data[prot_data == 0] <- NA
  #prot_data <- prot_data %>% rename(ID = accession)
  colnames(prot_data)[which(names(prot_data) == "accession")] <- "ID"
  colnames(prot_data)[which(names(prot_data) == "Protein Accession")] <- "ID"
  colnames(prot_data)[which(names(prot_data) == "Accession")] <- "ID"
  colnames(prot_data)[which(names(prot_data) == "DESCRIPTION")] <- "Description"
  colnames(prot_data)[which(names(prot_data) == "Description")] <- "Description"
  colnames(prot_data)[which(names(prot_data) == "description")] <- "Description" 
  colnames(prot_data)[which(names(prot_data) == "Protein Description")] <- "Description" 
  colnames(prot_data)[which(names(prot_data) == "Gene ID")] <- "genename"
  colnames(prot_data)[which(names(prot_data) == "gene_names")] <- "genename"
  colnames(prot_data)[which(names(prot_data) == "GeneName")] <- "genename"
  colnames(prot_data)[which(names(prot_data) == "GeneNames")] <- "genename"  
  colnames(prot_data)[which(names(prot_data) == "# AAs")] <- "protein_length_aa"
  colnames(prot_data)[which(names(prot_data) == "protein_length")] <- "protein_length_aa"
  prot_data %>% 
  write_xlsx(paste0(output_path,"tables/",name_df,"_selected_data.xlsx"))

# Remove Contaminant ####
# look for proteins labeled as contaminant and remove it

  text <- "contaminant"

  contaminants <- prot_data %>%
    dplyr::filter(stringr::str_detect(genename, text)) %>% 
    write_xlsx(paste0(output_path,"tables/",name_df,"_contaminants.xlsx"))


  prot_data <- prot_data %>%
    filter(!str_detect(genename, text))

  prot_data %>% 
    write_xlsx(paste0(output_path,"tables/",name_df,"_cleaned_data.xlsx"))
  
  # Update the dataframe in the global environment
  assign(name_df, prot_data, envir = .GlobalEnv)
return(prot_data)
}




'
# ENSEMBL column ####

# add translated symbol from zebrafish uniprot to zebrafish ensembl

values <- prot_data$accession

# Load biomaRt package
library(biomaRt)

# Connect to Ensembl database
ensembl <- useMart(biomart = "ensembl", dataset = "drerio_gene_ensembl", host = "https://www.ensembl.org")

# Retrieve Ensembl gene IDs for UniProt identifiers from zebrafish
zebrafish_mart <- useDataset("drerio_gene_ensembl", mart = ensembl)
zebrafish_genes <- biomaRt::getBM(attributes = c("ensembl_gene_id","uniprotswissprot", "external_gene_name", "zfin_id_symbol"), 
                         filters = "uniprotswissprot", 
                         values = values, 
                         mart = zebrafish_mart)

# clusterProfiler

ids <- bitr(prot_data$accession, fromType="UNIPROT", toType=c("ENSEMBL"), OrgDb="org.Dr.eg.db", drop = FALSE)
ids_non_ensembl <- ids %>% 
  filter(is.na(ENSEMBL))

ids_non_ensembl$UNIPROT %>% duplicated() %>% any()

# add column with converted gene names from zebrafish ensembl to human ensembl
'

