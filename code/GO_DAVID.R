# DAVID GO ####
#https://davidbioinformatics.nih.gov/content.jsp?file=DAVID_API.html

# libraries ####

source("../code/00_packages.R")
#source("../code/global_variables.R")

# Required packages
library(httr)
library(xml2)
library(dplyr)
library(purrr)


# DAVID API configuration
email <- "your.registered@email.com"  # Must be registered with DAVID
base_url <- "https://david.ncifcrf.gov/webservice/services/DAVIDWebService"


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
gene_list <- unlist(UPDOWN_dataframes_name, recursive = FALSE)#Repeat function twice for remove one level more of the list for GO


# GO functions ####

# DAVID GO API ####



#type: one of DAVID recognized gene types 
david_types <- c("AFFYMETRIX_3PRIME_IVT_ID", "AFFYMETRIX_EXON_GENE_ID","AFFYMETRIX_SNP_ID","AGILENT_CHIP_ID","AGILENT_ID","AGILENT_OLIGO_ID",
                 "ENSEMBL_GENE_ID", "ENSEMBL_TRANSCRIPT_ID","ENTREZ_GENE_ID","FLYBASE_GENE_ID","FLYBASE_TRANSCRIPT_ID","GENBANK_ACCESSION",
                 "GENPEPT_ACCESSION","GENOMIC_GI_ACCESSION","PROTEIN_GI_ACCESSION","ILLUMINA_ID","IPI_ID","MGI_ID","GENE_SYMBOL","PFAM_ID",
                 "PIR_ACCESSION","PIR_ID","PIR_NREF_ID","REFSEQ_GENOMIC","REFSEQ_MRNA","REFSEQ_PROTEIN","REFSEQ_RNA","RGD_ID","SGD_ID",
                 "TAIR_ID","UCSC_GENE_ID","UNIGENE","UNIPROT_ACCESSION","UNIPROT_ID","UNIREF100_ID","WORMBASE_GENE_ID","WORMPEP_ID","ZFIN_ID")



#annot, a list of desired annotation  categories separated by ","
annot_categories <- c(#Gene_Ontology
  "GOTERM_BP_1", "GOTERM_BP_2", "GOTERM_BP_3", "GOTERM_BP_4", "GOTERM_BP_5","GOTERM_BP_ALL", "GOTERM_BP_FAT",
  "GOTERM_CC_1", "GOTERM_CC_2", "GOTERM_CC_3", "GOTERM_CC_4", "GOTERM_CC_5","GOTERM_CC_ALL", "GOTERM_CC_FAT",
  "GOTERM_MF_1", "GOTERM_MF_2", "GOTERM_MF_3", "GOTERM_MF_4", "GOTERM_MF_5","GOTERM_MF_ALL", "GOTERM_MF_FAT",
  #Protein_Domains
  "BLOCKS_ID", "COG", "INTERPRO", "PDB_ID", "PFAM", "PIR_ALN", "PIR_HOMOLOGY_DOMAIN",
  "PIR_SUPERFAMILY", "PRINTS", "PRODOM", "PROSITE", "SCOP_ID", "SMART", "TIGRFAMS",
  #Pathways
  "BBID", "BIOCARTA", "EC_NUMBER", "KEGG_COMPOUND", "KEGG_PATHWAY", "KEGG_REACTION",
  #General Annotations
  "ALIAS_GENE_SYMBOL", "CHROMOSOME", "CYTOBAND", "GENE", "GENE_SYMBOL", "HOMOLOGOUS_GENE",
  "LL_SUMMARY", "OMIM_ID", "PIR_SUMMARY", "PROTEIN_MW", "REFSEQ_PRODUCT", "SEQUENCE_LENGTH",
  "SP_COMMENT",
  #Protein-Protein Interaction
  "CGAP_EST_QUARTILE", "CGAP_EST_RANK", "COG_ONTOLOGY", "PIR_SEQ_FEATURE", "SP_COMMENT_TYPE",
  "SP_PIR_KEYWORDS", "UP_SEQ_FEATURE",
  #Functional Categories
  "BIND", "DIP", "HIV_INTERACTION_CATEGORY", "HIV_INTERACTION", "MINT", "NCICB_CAPATHWAY",
  "TRANSFAC_ID",
  #Literature
  "GENERIF_SUMMARY", "HIV_INTERACTION_PUBMED_ID", "PUBMED_ID",
  #Disease
  "GENETIC_ASSOCIATION_DB_DISEASE", "OMIM_DISEASE")

#tool  = one of DAVID tool names
tool_type <- c("gene2gene", #Gene Functional Classification
               "term2term", #Funtional Annotation Clustering
               "summary", # Functional Annotation Summary
               "chartReport", #Functional Annotation Chart
               "annotationReport", #Functional Annotation Table
               "list", #Show Gene List Names in Batch
               "geneReport", #Gene Report
               "geneReportFull" #Gene Full Report
               )


type <- david_types[33]
annot <- annot_categories[c(6,12,18)]
tool <- tool_type[8]



# Función corregida para obtener resultados desde DAVID API
getDAVIDResults <- function(gene_list, user, idType = "OFFICIAL_GENE_SYMBOL", 
                            annotation = "GOTERM_BP_DIRECT", species = "Homo sapiens") {
  # URL de la API DAVID
  url <- "https://david.ncifcrf.gov/api.jsp"
  
  # Parámetros para la consulta POST incluyendo la clave de usuario DAVID
  params <- list(
    type = idType,
    ids = paste(gene_list, collapse = ","),
    tool = "chartReport",
    annot = annotation,
    user = user
  )
  
  # Realizar la solicitud POST
  response <- POST(url, body = params, encode = "form")
  
  # Verificar la respuesta del servidor
  if (status_code(response) == 200) {
    content_text <- content(response, "text", encoding = "UTF-8")
    if (nchar(content_text) == 0) {
      stop("Respuesta vacía desde DAVID. Revisa tu clave de usuario y parámetros.")
    }
    df <- read.delim(text = content_text, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    return(df)
  } else {
    stop("Error en la solicitud: ", status_code(response))
  }
}

# Ejemplo de uso con tu clave de usuario DAVID (user)
my_user_key <- "TU_CLAVE_DAVID"
results_list <- lapply(gene_list, getDAVIDResults, user = my_user_key)

# Mostrar resultados del primer vector
print(results_list[[1]])
