

# libraries ####
source("../code/00_packages.R")

# Load data ####

#prot_data <- readxl::read_xlsx(path = paste0(output_path,"tables/cleaned_data.xlsx"))

#Function takes as arguments 2 dataframes, one vector and one variable for the output path. One df with the protein data, other df with the experimental desing and a vector with the comparisons.

data_cleaning <- function(prot_data, Exp_design, comparisons, output_path, mnar_var = c("zero", "MinProb", "QRILC"), val_man = FALSE, value_imputation = 0.1){
  create_directories(paste0(output_path)) #Function from Globalvariables file to creates folder
  
  
  name_df <- deparse(substitute(prot_data))
# Make Unique ####
# Make unique names using the annotation in the "gene_names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
data_unique <- make_unique(prot_data, "genename","ID")
# Are there any duplicated names?
data_unique$Accession %>% duplicated() %>% any()


# Generate data_se element ####

# Generate a SummarizedExperiment object using an experimental design
LFQ_columns <- which(colnames(prot_data) %in% Exp_design$label)
# get LFQ column numbers
LFQ_columns

data_se <- make_se(data_unique, LFQ_columns, Exp_design)

# Generate a SummarizedExperiment object by parsing condition information from the LFQ column names
data_se_parsed <- make_se_parse(data_unique, LFQ_columns)

# Let's have a look at the SummarizedExperiment object
#data_se

# Filter on missing values ####
# The dataset contains proteins which are not quantified in all replicates. Some proteins are even only quantified in a single replicate.

# Plot a barplot of the protein identification overlap between samples
p<-plot_frequency(data_se)
print(p)
filename <- paste0(output_path,"figures/001_protein_identification_overlap_", name_df)# Name of output file. 
ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Adjust size according to your needs.
ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  

'This leaves our dataset with missing values, which need to be imputed.
However, this should not be done for proteins that contain too many
missing values. Therefore, we first filter out proteins that contain too
many missing values. This is done by setting the threshold for the
allowed number of missing values per condition in the filter_missvalfunction.'

# Less stringent filtering:
# Filter for proteins that are identified in 3 out of 4 replicates of at least one condition
# filters a proteomics dataset based on missing values. The dataset is filtered for  proteins that have a maximum of 'thr' missing values in at least one condition.
data_filt <- filter_missval(data_se, thr = 1)
save(data_filt, file = paste0(output_path,"RData/data_filt_", name_df,".RData"))

'After filtering, the number of identified proteins per sample can be
plotted as well as the overlap in identifications between samples.
'
p<- plot_numbers(data_filt)
print(p)
filename <- paste0(output_path,"figures/002_protein_per_sample_", name_df)# Name of output file. 
ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Adjust size according to your needs.
ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  

p<- plot_coverage(data_filt)
print(p)
filename <- paste0(output_path,"figures/003_protein_coverage_", name_df)# Name of output file. 
ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Adjust size according to your needs.
ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  

# Normalization ####
#The data is background corrected and normalized by variance stabilizing transformation (vsn).
# Normalize the data
data_norm <- normalize_vsn(data_filt)
save(data_norm, file = paste0(output_path,"RData/data_norm_", name_df,".RData"))

p <- summary(meanSdPlot(data_norm))
print(p)

p<- plot_normalization(data_filt, data_norm)
print(p)
filename <- paste0(output_path,"figures/004_normalized_data_", name_df)# Name of output file. 
ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Adjust size according to your needs.
ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  

# Imputation of missing values ####
' The remaining missing values in the dataset need to be imputed. The data
can be missing at random (MAR), for example if proteins are quantified
in some replicates but not in others. Data can also be missing not at
random (MNAR), for example if proteins are not quantified in specific
conditions (e.g. in the control samples). MNAR can indicate that
proteins are below the detection limit in specific samples, which could
be very well the case in proteomics experiments. For these different
conditions, different imputation methods have to be used, as described
in the MSnbase vignette and more specifically in the impute function
descriptions.

To explore the pattern of missing values in the data, a heatmap is
plotted indicating whether values are missing (0) or not (1). Only
proteins with at least one missing value are visualized. '

filename <- paste0(output_path,"figures/005_missing_values_", name_df)# Name of output file. 
# Step 1: Call the pdf command to start the plot
pdf(file = paste0(filename,pdf_extension),   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches

# Step 2: Create the plot with R code
p<-plot_missval(data_filt)
print(p)

# Step 3: Run dev.off() to create the file!
dev.off()

'This heatmap indicates that missing values are highly biased to specific
samples. To check whether missing values are biased to lower intense proteins, the densities and
cumulative fractions are plotted for proteins with and without missing
values.'

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
p<- plot_detect(data_filt)
print(p)
filename <- paste0(output_path,"figures/006_protein_imputation_", name_df)# Name of output file. 
ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Adjust size according to your needs.
ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  

'Indeed the proteins with missing values have on average low intensities.
This data (MNAR and close to the detection limit) should be imputed by a
left-censored imputation method, such as the quantile regression-based
left-censored function ("QRILC") or random draws from a left-shifted
distribution ("MinProb" and "man"). In contrast, MAR data should be
imputed with methods such as k-nearest neighbor ("knn") or maximum
likelihood ("MLE") functions. See the MSnbase vignette and more
specifically the impute function description for more information.'
set.seed(1)
# All possible imputation methods are printed in an error, if an invalid function name is given.
# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)

data_imp_QRILC <- impute(data_norm, fun = "QRILC")


## Advance Imputation method ####
### Mixed imputation on proteins (rows)
'(Taken from: https://bioconductor.org/packages/3.18/bioc/vignettes/DEP/inst/doc/MissingValues.html#roc-curves)
  One can also perform a mixed imputation on the proteins, which uses a MAR and MNAR imputation method on different subsets of proteins. 
  First, we have to define a logical vector defining the rows that are to be imputed with the MAR method. 
  Here, we consider a protein to have missing values not at random (MNAR) if it has missing values in all replicates of at least one condition.'
# Extract protein names with missing values 
# in all replicates of at least one condition




# Identificar proteínas MNAR
proteins_MNAR <- get_df_long(data_norm) %>%
  dplyr::mutate(across(where(is_character),as_factor)) %>% 
  group_by(name, condition) %>%
  summarize(MNAR_flag = sum(is.na(intensity)) / n() >= 0.5) %>%  # Cambiar criterio si es necesario
  filter(MNAR_flag) %>% 
  pull(name) %>% 
  unique()

'proteins_MNAR <- get_df_long(data_norm) %>%
  group_by(name, condition) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()


# 1. Create a lookup table that flags MNAR for each protein–condition group:
mnar_key <- get_df_long(data_norm) %>%
  mutate(across(where(is_character), as_factor)) %>%
  group_by(name, condition) %>%
  summarize(MNAR_flag = sum(is.na(intensity)) / n() >= 0.5, .groups = "drop") 

# 2. Pivot the lookup table into wide format:
mnar_key_wider <- mnar_key %>% 
  pivot_wider(names_from = condition, values_from = MNAR_flag) %>% 
  column_to_rownames(var = "name") %>% 
  as.matrix()
'

# Get a logical vector
MNAR <- rownames(data_norm) %in% proteins_MNAR

# Perform a mixed imputation
mixed_imputation <- impute(data_norm, 
                           fun = "mixed",
                           randna = MNAR, # we have to define MAR which is the opposite of MNAR
                           mar = "knn", # imputation function for MAR
                           #mnar = "zero") # imputation function for MNAR
                           mnar = mnar_var) # imputation function for MNAR

#val_man: Decide if you want to edit the mnar zero value.
if(mnar_var == "zero" && val_man == TRUE){
  assay(mixed_imputation)[assay(mixed_imputation) == 0] <- value_imputation
}
## Test for differential analysis for imputation methods
# We perform differential analysis on the different imputated data sets. The following datasets are compared: No imputation, knn imputation, MinProb imputation, Mixed imputation.

# Exporting results for data analysis ####

save(data_imp, file = paste0(output_path,"RData/data_imp_", name_df,".RData"))
save(data_imp_man, file = paste0(output_path,"RData/data_imp_man_", name_df,".RData"))
save(data_imp_knn, file = paste0(output_path,"RData/data_imp_knn_", name_df,".RData"))
save(mixed_imputation, file = paste0(output_path,"RData/mixed_imputation_", name_df,".RData"))
save(data_imp_QRILC, file = paste0(output_path,"RData/data_imp_QRILC_", name_df,".RData"))

load(paste0(output_path,"RData/data_imp_", name_df,".RData"))
load(paste0(output_path,"RData/data_imp_man_", name_df,".RData"))
load(paste0(output_path,"RData/data_imp_knn_", name_df,".RData"))
load(paste0(output_path,"RData/mixed_imputation_", name_df,".RData"))
load(paste0(output_path,"RData/data_norm_", name_df,".RData"))
load(paste0(output_path,"RData/data_imp_QRILC_", name_df,".RData"))

# Dynamically assign the data frame to the global environment
assign("data_imp", data_imp, envir = .GlobalEnv)
assign("data_imp_man", data_imp_man, envir = .GlobalEnv)
assign("data_imp_knn", data_imp_knn, envir = .GlobalEnv)
assign("mixed_imputation", mixed_imputation, envir = .GlobalEnv)
assign("data_imp_QRILC", mixed_imputation, envir = .GlobalEnv)



# Visualization of imputation effect
p<- plot_imputation(data_norm, data_imp, data_imp_man, data_imp_knn, mixed_imputation, data_imp_QRILC)
print(p)
filename <- paste0(output_path,"figures/007_protein_imputation_distribution_", name_df)# Name of output file. 
ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Adjust size according to your needs.
ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  

# For our dataset, knn and mixed imputation result in less identified differential expressed proteins compared to the no imputation and MinProb. No imputation results in the identification of the most differentially expressed proteins in our dataset with many proteins missing values.
# Note that the performance of the different imputation methods is data set-dependent. It is recommended to always carefully check the effect of filtering and data imputation on your results.

# PCA plots ####

## PCA plot ####
dep_analysis_norm <- analyze_dep(data_norm,type="manual",control=NULL, alpha=p_val,lfc=FC, test=comparisons ,design_formula=formula(~0+ condition))
dep_analysis_min <- analyze_dep(data_imp,type="manual",control=NULL, alpha=p_val,lfc=FC, test=comparisons ,design_formula=formula(~0+ condition))
dep_analysis_manual <- analyze_dep(data_imp_man,type="manual",control=NULL, alpha=p_val,lfc=FC, test=comparisons ,design_formula=formula(~0+ condition))
dep_analysis_knn <- analyze_dep(data_imp_knn,type="manual",control=NULL, alpha=p_val,lfc=FC, test=comparisons ,design_formula=formula(~0+ condition))
dep_analysis_mixed <- analyze_dep(mixed_imputation,type="manual",control=NULL, alpha=p_val,lfc=FC, test=comparisons ,design_formula=formula(~0+ condition))
dep_analysis_QRILC <- analyze_dep(data_imp_QRILC,type="manual",control=NULL, alpha=p_val,lfc=FC, test=comparisons ,design_formula=formula(~0+ condition))


#PCA MinProb
pca_minprob <- plot_pca(dep_analysis_min,x = 2, y = 1, n = (length(dep_analysis_min)), point_size = 4) + ggtitle("PCA MinProb",subtitle = paste0(length(dep_analysis_min), " variable proteins"))
print(pca_minprob)
filename <- paste0(output_path,"figures/008_PCA_MinProb_", name_df)# Name of output file. 
ggsave(paste0(filename,tiff_extension), pca_minprob, width = 8, height = 6, units = "in")  # Adjust size according to your needs.
ggsave(paste0(filename,pdf_extension), pca_minprob, width = 8, height = 6, units = "in") 
#PCA Manual
pca_manual <- plot_pca(dep_analysis_manual,x = 2, y = 1, n = (length(dep_analysis_manual)), point_size = 4)+ ggtitle("PCA Manual",subtitle = paste0(length(dep_analysis_manual), " variable proteins"))
print(pca_manual)
filename <- paste0(output_path,"figures/009_PCA_Manual_", name_df)# Name of output file. 
ggsave(paste0(filename,tiff_extension), pca_manual, width = 8, height = 6, units = "in")  # Adjust size according to your needs.
ggsave(paste0(filename,pdf_extension), pca_manual, width = 8, height = 6, units = "in") 
#PCA KNN
pca_knn <- plot_pca(dep_analysis_knn,x = 2, y = 1, n = (length(dep_analysis_knn)), point_size = 4)+ ggtitle("PCA KNN",subtitle = paste0(length(dep_analysis_knn), " variable proteins"))
print(pca_knn)
filename <- paste0(output_path,"figures/010_PCA_KNN_", name_df)# Name of output file. 
ggsave(paste0(filename,tiff_extension), pca_knn, width = 8, height = 6, units = "in")  # Adjust size according to your needs.
ggsave(paste0(filename,pdf_extension), pca_knn, width = 8, height = 6, units = "in") 
#PCA Mixed
pca_mixed <- plot_pca(dep_analysis_mixed,x = 2, y = 1, n = (length(dep_analysis_mixed)), point_size = 4)+ ggtitle("PCA Mixed Imputation",subtitle = paste0(length(dep_analysis_mixed), " variable proteins"))
print(pca_mixed)
filename <- paste0(output_path,"figures/011_PCA_Mixed_", name_df)# Name of output file. 
ggsave(paste0(filename,tiff_extension), pca_mixed, width = 8, height = 6, units = "in")  # Adjust size according to your needs.
ggsave(paste0(filename,pdf_extension), pca_mixed, width = 8, height = 6, units = "in") 
#PCA QRILC
pca_QRILC <- plot_pca(dep_analysis_QRILC,x = 2, y = 1, n = (length(dep_analysis_mixed)), point_size = 4)+ ggtitle("PCA QRILC",subtitle = paste0(length(dep_analysis_QRILC), " variable proteins"))
print(pca_QRILC)
filename <- paste0(output_path,"figures/012_PCA_QRILC_", name_df)# Name of output file. 
ggsave(paste0(filename,tiff_extension), pca_QRILC, width = 8, height = 6, units = "in")  # Adjust size according to your needs.
ggsave(paste0(filename,pdf_extension), pca_QRILC, width = 8, height = 6, units = "in") 



# PCA-plot Outlier identification ####
# Function to identify outliers within a group using robust PCA with error handling
run_pca_hubert_analysis <- function(pca_list, k = 5) {
  all_outliers <- data.frame(Sample = pca_list[[1]]$data$rowname)
  results <- list()
  
  for (pca_name in names(pca_list)) {
    pca_data <- pca_list[[pca_name]]$data
    
    robust_pca <- rrcov::PcaHubert(pca_data[, c("PC1", "PC2")], k = k)
    
    distances_df <- data.frame(
      Sample = pca_data$rowname,
      Condition = pca_data$condition,
      Score_Distance = robust_pca@sd,
      Orthogonal_Distance = robust_pca@od,
      Outlier = !robust_pca@flag
    )
    
    # Añadir al dataframe consolidado
    all_outliers <- all_outliers %>% 
      left_join(distances_df %>% select(Sample, Outlier), by = "Sample") %>% 
      rename(!!paste0("Outlier_", pca_name) := Outlier)
    
    # Generar gráfico
    plot_df <- distances_df
    p <- ggplot(plot_df, aes(x = Score_Distance, y = Orthogonal_Distance, label = Sample)) +
      geom_point(aes(color = Outlier), size = 3) +
      geom_text(hjust = 1.2, vjust = 0.5, size = 3) +
      scale_color_manual(values = c("black", "red")) +
      theme_minimal() +
      labs(title = paste("DD-plot PCAHubert:", pca_name),
           subtitle = "Outliers en rojo",
           x = "Score Distance",
           y = "Orthogonal Distance")
    print(p)
    filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_DD_Plot_", pca_name)
    ggsave(paste0(filename, ".tiff"), p, width = 8, height = 6, units = "in")
    ggsave(paste0(filename, ".pdf"), p, width = 8, height = 6, units = "in")
    
    assign("image_number", image_number + 1, envir = .GlobalEnv)
    
    results[[pca_name]] <- distances_df
  }
  
  # Guardar tabla consolidada de outliers
  write.xlsx(all_outliers, paste0(output_path, "tables/PCA_based_outliers.xlsx"))
  
  return(list("outlier_summary" = all_outliers, "detailed_results" = results))
}


# List of PCA objects
pca_objects <- list(
  pca_minprob = pca_minprob,
  pca_manual = pca_manual,
  pca_knn = pca_knn,
  pca_mixed = pca_mixed,
  pca_QRILC = pca_QRILC)
# Ejemplo de uso:
pca_results <- run_pca_hubert_analysis(pca_objects)

print(pca_results$outlier_summary)


save(pca_results, file = paste0(output_path,"RData/pca_outliers_results.RData"))


# Differential Expression ####
# Function that wraps around test_diff, add_rejections and get_results functions
DE_analysis <- function(se) {
  se %>% 
    test_diff(., type = "manual", test= comparisons) %>%
    add_rejections(., alpha = p_val, lfc = log2(FC)) %>% 
    get_results()
}

# DE analysis on no, knn, MinProb and mixed imputation
no_imputation_results <- DE_analysis(data_norm)
MinProb_imputation_results <- DE_analysis(data_imp)
Manual_imputation_results <- DE_analysis(data_imp_man)
knn_imputation_results <- DE_analysis(data_imp_knn)
mixed_imputation_results <- DE_analysis(mixed_imputation)
QRILC_imputation_results <- DE_analysis(data_imp_QRILC)


### Number of identified differentially expressed proteins
# As an initial parameter we look at the number of differentially expressed proteins identified (adjusted P ≤ 0.05 and Fold-Change > 1).

# Function to extract number of DE proteins
DE_prots <- function(results) {
  tibble(Dataset = gsub("_results", "", results),
         significant_proteins = get(results) %>% 
           filter(significant) %>% 
           nrow())
}

# Number of significant proteins
objects <- c("no_imputation_results", 
             "MinProb_imputation_results",
             "Manual_imputation_results",
             "knn_imputation_results",
             "mixed_imputation_results",
             "QRILC_imputation_results"
             )

dep_imputed_prots <- map_df(objects, DE_prots)
dep_imputed_prots$imputation_mnar <- mnar_var 
print(dep_imputed_prots)
dep_imputed_prots %>% 
  write_xlsx(paste0(output_path,"tables/dep_imputed_prots_", name_df,".xlsx"))

}
