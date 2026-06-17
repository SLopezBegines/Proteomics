# libraries ####

# source("../code/00_packages.R")
# Load data ####

# load(paste0(output_path,"RData/data_imp.RData"))
# load(paste0(output_path,"RData/data_imp_man.RData"))
# load(paste0(output_path,"RData/data_imp_knn.RData"))
# load(paste0(output_path,"RData/mixed_imputation.RData"))

# data_analysis() -------------------------------------------------------
# Differential enrichment analysis on the imputed SummarizedExperiment.
# Uses protein-wise linear models with empirical Bayes moderation (limma),
# wrapped by DEP::analyze_dep(). The model formula ~0 + condition fits a
# cell-means parameterisation so that manual contrasts map directly to
# biological comparisons without a reference level.
#
# Arguments:
#   imputation_file  SummarizedExperiment from data_cleaning() — typically
#                    mixed_splited_imputation
#   Exp_design       Experimental design table (condition, columns_to_rename …)
#   comparisons      Character vector of contrasts, e.g. "CTRL_vs_WT"
#   output_path      Base output directory
#
# Exports to global env: data_results, sig_adjusted_data, significative_data,
#   imputed_data_df, and one <comparison>_df per contrast.

data_analysis <- function(imputation_file, Exp_design, comparisons, output_path, adjusted_method = c("BH", "BY", "bonferroni", "holm", "fdr", "none")) {
  create_directories(output_path) # create tables/, figures/, RData/ subdirs

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


  # Fit limma model and test manual contrasts.
  # type = "manual"  → contrasts are taken from the 'test' vector verbatim
  # control = NULL   → no automatic control reference (handled by ~0 + condition)
  # alpha / lfc      → BH-adjusted p-value and log2FC thresholds from global_variables.R
  # design_formula   → cell-means model; each condition gets its own coefficient

  # Add Differential analysis if batch column is present in Exp_design, add batch as covariate in the design formula
  dep_analysis <- analyze_dep(imputation_file,
    type = "manual",
    control = NULL,
    alpha = p_val,
    lfc = FC,
    test = comparisons,
    design_formula = formula(~ 0 + condition)
  )


  imputed_data_df <- get_df_wide(imputation_file)
  # Generate a results table
  data_results <- get_results(dep_analysis)
  # Bind both df
  # label_names <- c("ID","Description","genename","protein_length_aa", "imputed","num_NAs")
  # label_names <- c("ID","Description","genename", "imputed","num_NAs")
  label_names <- c("ID", "Description", "genename", Exp_design$columns_to_rename)
  imputed_selected <- imputed_data_df[label_names]
  genenames <- imputed_selected$genename

  # Attempt to attach MNAR metadata (frac_NA and MNAR_flag per condition)
  # stored in the SummarizedExperiment by data_cleaning(). If the slot is
  # absent (e.g. when using a non-mixed imputation object), the merge is
  # skipped gracefully so downstream steps remain unaffected.
  imputed_selected <- tryCatch(
    {
      meta_df <- metadata(imputation_file)$proteins_MNAR %>%
        dplyr::filter(genename %in% genenames)
      merge(imputed_selected, meta_df, by = "genename", all = TRUE)
    },
    error = function(e) {
      warning("Skipping MNAR metadata merge (not available): ", e$message)
      imputed_selected
    }
  )

  data_results <- merge(data_results, imputed_selected, by = "ID", all = TRUE)

  ## add q-value corrections ####

  #  Supported `method` values in `p.adjust()`:
  #
  #  | Method | Type | When to use |
  #  |---|---|---|
  #  | `"BH"` | FDR (Benjamini-Hochberg) | Default; balanced power/control |
  #  | `"BY"` | FDR (Benjamini-Yekutieli) | When tests are positively dependent (correlated proteins) |
  #  | `"bonferroni"` | FWER | Very few tests; strict error control required |
  #  | `"holm"` | FWER | Step-down Bonferroni; slightly more powerful |
  #  | `"fdr"` | Alias for "BH" | — |
  #  | `"none"` | No correction | Only for visualization, never for publication |

  # Apply q-value per comparison
  for (comp in comparisons) {
    pval_col <- paste0(comp, "_p.val")
    qval_col <- paste0(comp, "_q.val")

    pvals <- data_results[[pval_col]]

    # qvalue() requires all p-values; NA causes failure
    valid <- !is.na(pvals)
    q_out <- qvalue(pvals[valid])

    data_results[[qval_col]] <- NA_real_
    data_results[[qval_col]][valid] <- q_out$qvalues

    adj_col <- paste0(comp, "_p.adj_holm") # example: Holm-Bonferroni
    data_results[[adj_col]] <- p.adjust(data_results[[pval_col]], method = adjusted_method)
  }


  # Add per-comparison regulation labels ####
  # For each contrast, add three columns:
  #   <comparison>_diffexpressed  "UP" | "DOWN" | "NO"  (nominal p-value threshold)
  #   <comparison>_dif_label      gene name if DE, NA otherwise (for volcano annotation)
  # Thresholds: ratio > FC & p.val < p_val → "UP"; ratio < -FC & p.val < p_val → "DOWN"
  for (i in 1:length(comparisons)) {
    #   p_val labels
    diff_pval_ <- paste(comparisons[i], "diffexpressed_pval", sep = "_")
    ratio <- paste(comparisons[i], "ratio", sep = "_")
    pval <- paste(comparisons[i], "p.val", sep = "_")
    label_pval <- paste(comparisons[i], "dif_label_pval", sep = "_")

    # add a column of NAs
    data_results[diff_pval_] <- "NO"
    # if log2Foldchange > FC and pvalue < p_val, set as "UP"
    data_results[diff_pval_][data_results[ratio] > FC & data_results[pval] < p_val] <- "UP"
    # if log2Foldchange < -FC and pvalue < p_val, set as "DOWN"
    data_results[diff_pval_][data_results[ratio] < -FC & data_results[pval] < p_val] <- "DOWN"
    # Now write down the name of genes beside the points...
    # Create a new column "dif_label_dataset" to data_results, that will contain the name of genes differentially expressed (NA in case they are not)
    data_results[label_pval] <- NA
    # add a column for names of differentialy expressed genes to label points in vulcano plot
    data_results[label_pval][data_results[diff_pval_] != "NO"] <- data_results$name[data_results[diff_pval_] != "NO"]


    # q_value labels
    diff_qval_ <- paste(comparisons[i], "diffexpressed_qval", sep = "_")
    qval <- paste(comparisons[i], "q.val", sep = "_")
    label_qval <- paste(comparisons[i], "dif_label_qval", sep = "_")


    # add a column of NAs
    data_results[diff_qval_] <- "NO"
    # if log2Foldchange > FC and pvalue < p_val, set as "UP"
    data_results[diff_qval_][data_results[ratio] > FC & data_results[qval] < p_val] <- "UP"
    # if log2Foldchange < -FC and pvalue < p_val, set as "DOWN"
    data_results[diff_qval_][data_results[ratio] < -FC & data_results[qval] < p_val] <- "DOWN"
    # Now write down the name of genes beside the points...

    # Create a new column "dif_label_dataset" to data_results, that will contain the name of genes differentially expressed (NA in case they are not)
    data_results[label_qval] <- NA
    # add a column for names of differentialy expressed genes to label points in vulcano plot
    data_results[label_qval][data_results[diff_qval_] != "NO"] <- data_results$name[data_results[diff_qval_] != "NO"]


    # padj_holm_value labels
    diff_holmval_ <- paste(comparisons[i], "diffexpressed_padj_holm_pval", sep = "_")
    holmval <- paste(comparisons[i], "p.adj_holm", sep = "_")
    label_holmval <- paste(comparisons[i], "dif_label_holmval", sep = "_")


    # add a column of NAs
    data_results[diff_holmval_] <- "NO"
    # if log2Foldchange > FC and pvalue < p_val, set as "UP"
    data_results[diff_holmval_][data_results[ratio] > FC & data_results[holmval] < p_val] <- "UP"
    # if log2Foldchange < -FC and pvalue < p_val, set as "DOWN"
    data_results[diff_holmval_][data_results[ratio] < -FC & data_results[holmval] < p_val] <- "DOWN"
    # Now write down the name of genes beside the points...

    # Create a new column "dif_label_dataset" to data_results, that will contain the name of genes differentially expressed (NA in case they are not)
    data_results[label_holmval] <- NA
    # add a column for names of differentialy expressed genes to label points in vulcano plot
    data_results[label_holmval][data_results[diff_holmval_] != "NO"] <- data_results$name[data_results[diff_holmval_] != "NO"]
  }

  # Add suffix "lfq_" to specified column names
  data_results <- data_results %>%
    rename_with(~ paste0("lfq_", .x), any_of(Exp_design$columns_to_rename))

  # Get col names
  column_names <- names(data_results)

  # Create a new null column called significance. This will be TRUE if in at least one comparison FC and pval pass threshold.
  data_results$significance_pval <- NA
  data_results$significance_qval <- NA
  data_results$significance_holmval <- NA

  # After the loop, set significance to TRUE if at least one diff_ column contains "UP" or "DOWN"
  data_results <- data_results %>%
    rowwise() %>%
    mutate(
      significance_pval = any(c_across(ends_with("_diffexpressed_pval")) %in% c("UP", "DOWN")),
      significance_qval = any(c_across(ends_with("_diffexpressed_qval")) %in% c("UP", "DOWN")),
      significance_holmval = any(c_across(ends_with("_diffexpressed_padj_holm_pval")) %in% c("UP", "DOWN")),
    ) %>%
    ungroup()

  # Composite score = log2FC × −log10(p-value).
  # Combines effect size and significance into a single ranking metric,
  # equivalent to the signed volcano distance from the origin.
  for (i in comparisons) {
    ratio <- paste0(i, "_ratio")
    p.val <- paste0(i, "_p.val")
    score <- paste0(i, "_Score")
    data_results <- data_results %>%
      mutate(!!score := !!sym(ratio) * -log10(!!sym(p.val)))
  }

  # Get col names
  column_names <- names(data_results)
  # Sort cols by termination
  # sorted_column_names <- c("name", "ID", "Description", "protein_length_aa", "significant",
  sorted_column_names <- c(
    "name", "ID", "Description", "significant", "significance_pval", "significance_qval", "significance_holmval",
    column_names[grep("_ratio$", column_names)],
    column_names[grep("_p.val$", column_names)],
    column_names[grep("_p.adj$", column_names)],
    column_names[grep("_q.val$", column_names)],
    column_names[grep("_p.adj_holm$", column_names)],
    column_names[grep("_Score$", column_names)],
    column_names[grep("_centered$", column_names)],
    column_names[grep("_diffexpressed_pval$", column_names)],
    column_names[grep("_dif_label_pval$", column_names)],
    column_names[grep("_diffexpressed_qval$", column_names)],
    column_names[grep("_dif_label_qval$", column_names)],
    column_names[grep("_diffexpressed_padj_holm_pval$", column_names)],
    column_names[grep("dif_label_holmval$", column_names)],
    column_names[grep("lfq_", column_names)]
  )

  # Agregar columnas restantes que no fueron especificadas
  remaining_cols <- setdiff(column_names, sorted_column_names)
  sorted_column_names <- c(sorted_column_names, remaining_cols)


  data_results <- data_results[, sorted_column_names]


  # Get list of differential gene by comparison
  differential_list <- list()
  for (i in comparisons) {
    diffgenes <- paste0(i, "_diffexpressed_pval")
    ratio <- paste0(i, "_ratio")
    p.val <- paste0(i, "_p.val")
    adj.pval <- paste0(i, "_p.adj")
    q.val <- paste0(i, "_q.val")
    holm.val <- paste0(i, "_p.adj_holm")
    score <- paste0(i, "_Score")


    diff <- data_results %>%
      dplyr::filter(.data[[diffgenes]] == "UP" | .data[[diffgenes]] == "DOWN") %>%
      dplyr::select(name, ID, Description, all_of(diffgenes), all_of(ratio), all_of(p.val), all_of(adj.pval), all_of(q.val), all_of(score)) %>%
      dplyr::arrange(desc(!!sym(score)))

    differential_list[[i]] <- diff # Almacenar los resultados en la lista

    # Crear un nuevo data frame con un nombre basado en la variable original
    new_df_name <- paste(i, "df", sep = "_")

    assign(new_df_name, as.data.frame(diff), envir = .GlobalEnv)

    diff %>%
      write_xlsx(paste0(output_path, "tables/", new_df_name, ".xlsx"))
  }

  ## p-value distribution plot ####
  pval_cols <- grep("_p.val$", names(data_results), value = TRUE)
  pval_data <- data_results %>%
    select(all_of(pval_cols)) %>%
    pivot_longer(cols = everything(), names_to = "comparison", values_to = "p_value") %>%
    mutate(comparison = str_replace(comparison, "_p.val$", ""))


  pval_plot <- ggplot(pval_data, aes(x = p_value)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "black") +
    facet_wrap(~comparison, scales = "free_y") +
    theme_minimal() +
    labs(title = "P-value Distribution", x = "P-value", y = "Frequency")


  ## p.adj distribution plot ####
  p.adj_cols <- grep("_p.adj$", names(data_results), value = TRUE)
  p.adj_data <- data_results %>%
    select(all_of(p.adj_cols)) %>%
    pivot_longer(cols = everything(), names_to = "comparison", values_to = "p.adj") %>%
    mutate(comparison = str_replace(comparison, "_p.adj$", ""))


  p.adj_plot <- ggplot(p.adj_data, aes(x = p.adj)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "black") +
    facet_wrap(~comparison, scales = "free_y") +
    theme_minimal() +
    labs(title = "Adjusted p-value Distribution", x = "Adjusted p-value", y = "Frequency")

  ## q-value distribution plot ####
  qval_cols <- grep("_q.val$", names(data_results), value = TRUE)
  qval_data <- data_results %>%
    select(all_of(qval_cols)) %>%
    pivot_longer(cols = everything(), names_to = "comparison", values_to = "q_value") %>%
    mutate(comparison = str_replace(comparison, "_q.val$", ""))


  qval_plot <- ggplot(qval_data, aes(x = q_value)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "black") +
    facet_wrap(~comparison, scales = "free_y") +
    theme_minimal() +
    labs(title = "Q-value Distribution", x = "Q-value", y = "Frequency")

  ## holm adjusted-pvalue distribution plot ####
  holm_cols <- grep("_p.adj_holm$", names(data_results), value = TRUE)
  holm_data <- data_results %>%
    select(all_of(holm_cols)) %>%
    pivot_longer(cols = everything(), names_to = "comparison", values_to = "holm_adjusted_pvalue") %>%
    mutate(comparison = str_replace(comparison, "_p.adj_holm$", ""))


  holm_plot <- ggplot(holm_data, aes(x = holm_adjusted_pvalue)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "black") +
    facet_wrap(~comparison, scales = "free_y") +
    theme_minimal() +
    labs(title = paste0("Adjusted ", adjusted_method, " method p-value Distribution"), x = paste0("Adjusted ", adjusted_method, " p-value"), y = "Frequency")


  adjusted_plot <- pval_plot / p.adj_plot / qval_plot / holm_plot

  save_plot("pq_value_distibution", adjusted_plot,
    output_dir = paste0(output_path, "figures")
  )

  # Number of significant proteins by pval/adj.pval/qval ####
  n_significant_adjp <- data_results %>%
    filter(significant == TRUE) %>%
    nrow()
  n_significant_pval <- data_results %>%
    filter(significance_pval == TRUE) %>%
    nrow()
  n_significant_qval <- data_results %>%
    filter(significance_qval == TRUE) %>%
    nrow()
  n_significant_holmval <- data_results %>%
    filter(significance_holmval == TRUE) %>%
    nrow()

  # Dataframe for number of significant proteins according pval/adj.pval/qval
  n_significant <- data.frame(
    c("p-value", "adjusted p-value", "q-value", paste0(adjusted_method, " adjusted p-value")),
    c(n_significant_pval, n_significant_adjp, n_significant_qval, n_significant_holmval)
  )
  colnames(n_significant) <- c("N significant proteins by", "Number of proteins")
  n_significant |>    write_xlsx(paste0(output_path, "tables/number_significant_proteins_by p-value.xlsx"))
  print(n_significant)

# Export results ####
  significant_adjusted_data <- data_results %>% filter(significant)
  significant_pval_data <- data_results %>% filter(significance_pval)
  significant_qval_data <- data_results %>% filter(significance_qval)
  significant_holmval_data <- data_results %>% filter(significance_holmval)

  save(dep_analysis, file = paste0(output_path, "RData/dep_analysis.RData"))
  # All data
  data_results %>%
    write_xlsx(paste0(output_path, "tables/data_results.xlsx"))
  # Significant by p.adj
  significant_adjusted_data %>%
    write_xlsx(paste0(output_path, "tables/significant_adjusted_data.xlsx"))
  # Significatn by p.val
  significant_pval_data %>%
    write_xlsx(paste0(output_path, "tables/significant_pval_data.xlsx"))
  # Significant by q.val
  significant_qval_data %>%
    write_xlsx(paste0(output_path, "tables/significant_qval_data.xlsx"))
  # Significant by holm p.val
  significant_holmval_data %>%
    write_xlsx(paste0(output_path, "tables/significant_holmval_data.xlsx"))

  # Dynamically assign the data frame to the global environment
  assign("data_results", data_results, envir = .GlobalEnv)
  assign("significant_adjusted_data", significant_adjusted_data, envir = .GlobalEnv)
  assign("significant_pval_data", significant_pval_data, envir = .GlobalEnv)
  assign("significant_qval_data", significant_qval_data, envir = .GlobalEnv)
  assign("significant_holmval_data", significant_holmval_data, envir = .GlobalEnv)
  assign("imputed_data_df", imputed_data_df, envir = .GlobalEnv)
}
