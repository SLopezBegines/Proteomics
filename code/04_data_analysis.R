

# libraries ####

source("../code/00_packages.R")
# Load data ####

#load(paste0(output_path,"RData/data_imp.RData"))
# load(paste0(output_path,"RData/data_imp_man.RData"))
# load(paste0(output_path,"RData/data_imp_knn.RData"))
# load(paste0(output_path,"RData/mixed_imputation.RData"))

data_analysis <- function(imputation_file, Exp_design, comparisons, output_path){
  create_directories(output_path) #Function from Globalvariables file to creates folder
  
  
  name_df <- deparse(substitute(imputation_file))
# Differential Expression Analysis ####

' Protein-wise linear models combined with empirical Bayes statistics are
used for the differential enrichment analysis (or differential
expression analysis). The test_diff function introduced here uses limma
and automatically generates the contrasts to be tested. For the
contrasts generation, the control sample has to be specified.
Additionally, the types of contrasts to be produced need to be
indicated, allowing the generation of all possible comparisons ("all")
or the generation of contrasts of every sample versus control
("control"). Alternatively, the user can manually specify the contrasts
to be tested (type = "manual"), which need to be specified in the
argument test.'

  # Differential enrichment analysis  based on linear models and empherical Bayes statistics
  
  # Test every sample versus control
    dep_analysis <- analyze_dep(imputation_file,
                                type="manual",
                                control=NULL, 
                                alpha=p_val,
                                lfc=FC, 
                                test=comparisons ,
                                design_formula=formula(~0+ condition))
    
  
  imputed_data_df <- get_df_wide(imputation_file)
  # Generate a results table
  data_results <- get_results(dep_analysis)
  # Bind both df
  #label_names <- c("ID","Description","genename","protein_length_aa", "imputed","num_NAs")
  #label_names <- c("ID","Description","genename", "imputed","num_NAs")
  label_names <- c("ID","Description","genename", Exp_design$columns_to_rename)
  imputed_selected <- imputed_data_df[label_names]
  genenames <- imputed_selected$genename
  
  imputed_selected <- tryCatch({
    
    # 1. pull out the MNAR‐metadata (may error if it's not there)
    meta_df <- metadata(imputation_file)$proteins_MNAR %>% 
      dplyr::filter(genename %in% genenames)
    
    # 2. merge it in
    merge(imputed_selected, meta_df, by = "genename", all = TRUE)
    
  }, error = function(e) {
    
    # if anything went wrong, warn and return unaltered
    warning("Skipping metadata merge (not available): ", e$message)
    imputed_selected
    
  })
  
  data_results <- merge(data_results, imputed_selected, by = "ID", all = TRUE)
  
  # Number of significant proteins
  data_results %>% filter(significant) %>% nrow()
  #data_results
  
  sig_adjusted_data <- data_results%>% filter(significant)
  

  
  # Add significance by comparison ####

  # Create one column for each comparison were it is indicated if the protein is Up-regulated or Down-regulated. 
  #Setting up dataframe
  for (i in 1:length(comparisons)){
    diff_ <- paste(comparisons[i],"diffexpressed", sep="_")
    ratio <-paste(comparisons[i],"ratio", sep="_")
    pval <- paste(comparisons[i],"p.val", sep="_")
    label <- paste(comparisons[i],"dif_label",sep="_") 
    # add a column of NAs
    data_results[diff_] <- "NO"
    # if log2Foldchange > FC and pvalue < p_val, set as "UP" 
    data_results[diff_][data_results[ratio] > FC & data_results[pval] < p_val] <- "UP"
    # if log2Foldchange < -FC and pvalue < p_val, set as "DOWN"
    data_results[diff_][data_results[ratio] < -FC & data_results[pval] < p_val] <- "DOWN"
    # Now write down the name of genes beside the points...
    # Create a new column "dif_label_dataset" to data_results, that will contain the name of genes differentially expressed (NA in case they are not)
    data_results[label] <- NA
    #add a column for names of differentialy expressed genes to label points in vulcano plot
    data_results[label][data_results[diff_] != "NO"] <- data_results$name[data_results[diff_] != "NO"]
    }
  
  # Add suffix "lfq_" to specified column names
  data_results <- data_results %>%
    rename_with(~ paste0("lfq_", .x), any_of(Exp_design$columns_to_rename))
  
  # Get col names
  column_names <- names(data_results)
  
  # Sort cols by termination
  #sorted_column_names <- c("name", "ID", "Description", "protein_length_aa", "significant",
  sorted_column_names <- c("name", "ID", "Description", "significant",
                           column_names[grep("_ratio$", column_names)],
                           column_names[grep("_p.val$", column_names)],
                           column_names[grep("_p.adj$", column_names)],
                           column_names[grep("_centered$", column_names)],
                           column_names[grep("_diffexpressed$", column_names)],
                           column_names[grep("_dif_label$", column_names)],
                           column_names[grep("_significant", column_names)],
                           column_names[grep("lfq_", column_names)])
  
  # Agregar columnas restantes que no fueron especificadas
  remaining_cols <- setdiff(column_names, sorted_column_names)
  sorted_column_names <- c(sorted_column_names, remaining_cols)
  
  
  data_results <- data_results[, sorted_column_names]
 

  # Create a new null column called significance. This will be TRUE if in at least one comparison FC and pval pass threshold.
  data_results$significance <- NA
  
  # After the loop, set significance to TRUE if at least one diff_ column contains "UP" or "DOWN"
  data_results <- data_results %>%
    rowwise() %>%
    mutate(significance = any(c_across(ends_with("_diffexpressed")) %in% c("UP", "DOWN"))) %>%
    ungroup()
  
  for (i in comparisons){
    ratio <- paste0(i,"_ratio")
    p.val <- paste0(i,"_p.val")
    score <- paste0(i, "_Score")
    data_results <- data_results %>% 
      mutate(!!score := !!sym(ratio) * -log10(!!sym(p.val)))
  }
  
  # Get col names
  column_names <- names(data_results)
  # Sort cols by termination
  #sorted_column_names <- c("name", "ID", "Description", "protein_length_aa", "significant","significance", 
  sorted_column_names <- c("name", "ID", "Description", "significant","significance",                         
                           column_names[grep("_ratio$", column_names)],
                           column_names[grep("_p.val$", column_names)],
                           column_names[grep("_p.adj$", column_names)],
                           column_names[grep("_Score$", column_names)],
                           column_names[grep("_centered$", column_names)],
                           column_names[grep("_diffexpressed$", column_names)],
                           column_names[grep("_dif_label$", column_names)],
                           #column_names[grep("_significant", column_names)],
                           column_names[grep("lfq_", column_names)])
  
  # Agregar columnas restantes que no fueron especificadas
  remaining_cols <- setdiff(column_names, sorted_column_names)
  sorted_column_names <- c(sorted_column_names, remaining_cols)
  
  data_results <- data_results[, sorted_column_names]
  


  #Get list of differential gene by comparison
  differential_list <- list()
  for (i in comparisons){
    diffgenes <- paste0(i,"_diffexpressed")
    ratio <- paste0(i,"_ratio")
    p.val <- paste0(i,"_p.val")
    adj.pval <- paste0(i,"_p.adj")
    score <- paste0(i,"_Score")
    
    
    diff <- data_results %>% 
      dplyr::filter(.data[[diffgenes]]=="UP" | .data[[diffgenes]]== "DOWN") %>% 
      dplyr::select(name, ID,Description, all_of(diffgenes), all_of(ratio), all_of(p.val), all_of(adj.pval), all_of(score)) %>% 
      dplyr::arrange(desc(!!sym(score)))
    
    differential_list[[i]] <- diff  # Almacenar los resultados en la lista
    
    # Crear un nuevo data frame con un nombre basado en la variable original
    new_df_name <- paste(i, "df", sep = "_")
    
    assign(new_df_name, as.data.frame(diff), envir = .GlobalEnv)
    
    diff %>% 
      write_xlsx(paste0(output_path,"tables/", new_df_name, ".xlsx"))
  }
  
  data_results %>% filter(significance == TRUE) %>% summarise("Total differential" = n())
  
  significative_data <- data_results %>% filter(significance == TRUE)

  
  sig_adjusted_data <- data_results%>% filter(significant == TRUE)
  
  # Export results ####
  
  save(dep_analysis, file = paste0(output_path,"RData/dep_analysis.RData"))
  
  data_results %>% 
    write_xlsx(paste0(output_path,"tables/data_results.xlsx"))
  sig_adjusted_data %>% 
    write_xlsx(paste0(output_path,"tables/sig_adjusted_data.xlsx"))
  significative_data %>% 
    write_xlsx(paste0(output_path,"tables/significative_data.xlsx"))

  
  # Dynamically assign the data frame to the global environment
  assign("data_results", data_results, envir = .GlobalEnv)
  assign("sig_adjusted_data", sig_adjusted_data, envir = .GlobalEnv)
  assign("significative_data", significative_data, envir = .GlobalEnv)
  assign("imputed_data_df", imputed_data_df, envir = .GlobalEnv)
  
  
  
  
}
