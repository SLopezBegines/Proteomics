# libraries ####
source("../code/00_packages.R")

# Load data ####

# prot_data <- readxl::read_xlsx(path = paste0(output_path,"tables/cleaned_data.xlsx"))

# Function takes as arguments 2 dataframes, one vector and one variable for the output path. One df with the protein data, other df with the experimental desing and a vector with the comparisons.

data_cleaning <- function(prot_data, # input dataset
                          Exp_design, # Experimental matrix
                          comparisons, # Comparisons to be made
                          output_path, # output path
                          fraction_NA = c(0.4, 0.5, 0.6),
                          factor_SD_impute = c(0.05, 0.1, 0.2),
                          mnar_var = c("zero", "MinProb", "QRILC")) {
  create_directories(paste0(output_path)) # Function from Globalvariables file to creates folder


  name_df <- deparse(substitute(prot_data))
  # Make Unique ####
  # Make unique names using the annotation in the "gene_names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
  data_unique <- make_unique(prot_data, "genename", "ID")
  # Are there any duplicated names?
  data_unique$Accession %>%
    duplicated() %>%
    any()


  # Generate data_se element ####

  # Generate a SummarizedExperiment object using an experimental design
  LFQ_columns <- which(colnames(prot_data) %in% Exp_design$label)
  # get LFQ column numbers
  LFQ_columns

  data_se <- make_se(data_unique, LFQ_columns, Exp_design)

  # Generate a SummarizedExperiment object by parsing condition information from the LFQ column names
  data_se_parsed <- make_se_parse(data_unique, LFQ_columns)

  # Let's have a look at the SummarizedExperiment object
  # data_se

  # Filter on missing values ####
  # The dataset contains proteins which are not quantified in all replicates. Some proteins are even only quantified in a single replicate.

  # Plot a barplot of the protein identification overlap between samples
  p1 <- plot_frequency(data_se)
  print(p1)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_protein_identification_overlap_", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), p1, width = 8, height = 6, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), p1, width = 8, height = 6, units = "in")
  image_number <<- image_number + 1
  "This leaves our dataset with missing values, which need to be imputed.
However, this should not be done for proteins that contain too many
missing values. Therefore, we first filter out proteins that contain too
many missing values. This is done by setting the threshold for the
allowed number of missing values per condition in the filter_missvalfunction."

  # Less stringent filtering:
  # Filter for proteins that are identified in 3 out of 4 replicates of at least one condition
  # filters a proteomics dataset based on missing values. The dataset is filtered for  proteins that have a maximum of 'thr' missing values in at least one condition.
  data_filt <- filter_missval(data_se, thr = 1)
  save(data_filt, file = paste0(output_path, "RData/data_filt_", name_df, ".RData"))

  "After filtering, the number of identified proteins per sample can be
plotted as well as the overlap in identifications between samples.
"
  p2 <- plot_numbers(data_filt)
  print(p2)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_protein_per_sample_", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), p2, width = 8, height = 6, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), p2, width = 8, height = 6, units = "in")
  image_number <<- image_number + 1
  p3 <- plot_coverage(data_filt)
  print(p3)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_protein_coverage_", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), p3, width = 8, height = 6, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), p3, width = 8, height = 6, units = "in")
  image_number <<- image_number + 1
  # Normalization ####
  # The data is background corrected and normalized by variance stabilizing transformation (vsn).
  # Normalize the data
  data_norm <- normalize_vsn(data_filt)
  save(data_norm, file = paste0(output_path, "RData/data_norm_", name_df, ".RData"))

  p <- summary(DEP::meanSdPlot(data_norm))
  print(p)

  p4 <- plot_normalization(data_filt, data_norm)
  print(p4)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_normalized_data_", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), p4, width = 8, height = 6, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), p4, width = 8, height = 6, units = "in")
  image_number <<- image_number + 1


  p <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2) + plot_annotation(
    title = "Data overview",
    subtitle = paste0("Dataset: ", name_df, " | Proteins: ", nrow(data_filt), " | Samples: ", ncol(data_filt)),
    caption = "Figures generated with the DEP package",
    tag_levels = "A"
  )
  print(p)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_QC_data_overview_", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), p, width = 16, height = 12, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), p, width = 16, height = 12, units = "in")
  image_number <<- image_number + 1

  "The normalization can be assessed by plotting the standard deviation
versus the mean intensity for each protein before and after normalization.
Additionally, a PCA plot can be used to assess whether the replicates
cluster together and whether there are any outliers. Both plots should
be generated before and after normalization to assess the effect of the
normalization."

  # Visual diagnostics of normalization ####
  ## Plot PCA and meanSdPlot before and after normalization ####
  plot_vsn_diagnostics <- function(se_raw, se_norm,
                                   group_col = "condition",
                                   output_path = output_path,
                                   file_prefix = "vsn_diagnostics") {
    # Remove "/" in output_path

    output_path_short <- gsub("/$", "", output_path)
    # Extraer matrices de expresión
    mat_raw <- assay(se_raw)
    mat_norm <- assay(se_norm)

    # Eliminar filas con NA para PCA (pero no para meanSdPlot)
    mat_raw_pca <- mat_raw[complete.cases(mat_raw), ]
    mat_norm_pca <- mat_norm[complete.cases(mat_norm), ]

    # Extraer metadatos
    meta <- as.data.frame(colData(se_raw))

    # PCA antes de normalización
    pca_raw <- prcomp(t(mat_raw_pca), scale. = TRUE)
    pca_raw_df <- as.data.frame(pca_raw$x[, 1:2])
    pca_raw_df$condition <- meta$condition
    pca_raw_df$replicate <- as.factor(meta$replicate)

    pca_raw_plot <- ggplot(pca_raw_df, aes(x = PC1, y = PC2, color = condition, shape = replicate)) +
      geom_point(size = 4, alpha = 0.8) +
      labs(title = "PCA - Before VSN") +
      theme_minimal() +
      theme(legend.position = "right")

    # PCA después de normalización
    pca_norm <- prcomp(t(mat_norm_pca), scale. = TRUE)
    pca_norm_df <- as.data.frame(pca_norm$x[, 1:2])
    pca_norm_df$condition <- meta$condition
    pca_norm_df$replicate <- as.factor(meta$replicate)

    pca_norm_plot <- ggplot(pca_norm_df, aes(x = PC1, y = PC2, color = condition, shape = replicate)) +
      geom_point(size = 4, alpha = 0.8) +
      labs(title = "PCA - After VSN") +
      theme_minimal() +
      theme(legend.position = "right")

    # meanSdPlot antes (usando px/py)
    mean_sd_raw_obj <- DEP::meanSdPlot(mat_raw, plot = FALSE)
    mean_sd_raw_data <- mean_sd_raw_obj$gg$data
    colnames(mean_sd_raw_data) <- c("x", "y") # px → x, py → y

    mean_sd_raw_plot <- ggplot(mean_sd_raw_data, aes(x = x, y = y)) +
      geom_point(alpha = 0.7) +
      geom_smooth(se = FALSE, color = "blue") +
      labs(title = "meanSdPlot - Before VSN", x = "rank(mean)", y = "SD") +
      theme_minimal()

    # meanSdPlot después
    mean_sd_norm_obj <- DEP::meanSdPlot(mat_norm, plot = FALSE)
    mean_sd_norm_data <- mean_sd_norm_obj$gg$data
    colnames(mean_sd_norm_data) <- c("x", "y")

    mean_sd_norm_plot <- ggplot(mean_sd_norm_data, aes(x = x, y = y)) +
      geom_point(alpha = 0.7) +
      geom_smooth(se = FALSE, color = "blue") +
      labs(title = "meanSdPlot - After VSN", x = "rank(mean)", y = "SD") +
      theme_minimal()

    # Componer figura
    final_plot <- (pca_raw_plot | pca_norm_plot) /
      (mean_sd_raw_plot | mean_sd_norm_plot)
    print(final_plot)
    # Guardar si se especifica output_path
    if (!is.null(output_path)) {
      pdf_file <- file.path(output_path_short, "figures", paste0(sprintf("%02d", image_number), "_", file_prefix, ".pdf"))
      tiff_file <- file.path(output_path_short, "figures", paste0(sprintf("%02d", image_number), "_", file_prefix, ".tiff"))

      ggsave(pdf_file, final_plot, width = 12, height = 10)
      ggsave(tiff_file, final_plot, width = 12, height = 10, dpi = 300)
      message("✅ Graph stored in ", output_path)
    }

    return(final_plot)
  }


  # Asumiendo que ya tienes 'se_raw' y 'se_norm'
  # se_raw <- filt
  # se_norm <- normalize_vsn(filt)
  plot_vsn_diagnostics(
    se_raw = data_filt,
    se_norm = data_norm,
    output_path = output_path,
    file_prefix = "Normalization_diagnosis"
  )
  image_number <<- image_number + 1
  # Imputation of missing values ####
  " The remaining missing values in the dataset need to be imputed. The data
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
proteins with at least one missing value are visualized. "

  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_missing_values_", name_df) # Name of output file.
  # Step 1: Call the pdf command to start the plot
  pdf(
    file = paste0(filename, pdf_extension), # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4
  ) # The height of the plot in inches

  # Step 2: Create the plot with R code
  p <- plot_missval(data_filt)
  print(p)

  # Step 3: Run dev.off() to create the file!
  dev.off()
  print(p)
  image_number <<- image_number + 1
  "This heatmap indicates that missing values are highly biased to specific
samples. To check whether missing values are biased to lower intense proteins, the densities and
cumulative fractions are plotted for proteins with and without missing
values."

  # Plot intensity distributions and cumulative fraction of proteins with and without missing values
  p <- plot_detect(data_filt)
  print(p)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_protein_imputation_", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), p, width = 8, height = 6, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), p, width = 8, height = 6, units = "in")
  image_number <<- image_number + 1
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
  data_imp_man_gauss <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

  # Impute missing data using the k-nearest neighbour approach (for MAR)
  data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)

  data_imp_QRILC <- impute(data_norm, fun = "QRILC")

  # Imputation of manual value
  # Global minimum for left-censored distribution
  global_min <- min(assay(data_norm), na.rm = TRUE)
  # Impute missing values with a value that is 1 SD below the global minimum
  value_impute <- global_min - sd(assay(data_norm), na.rm = TRUE) * factor_SD_impute
  # Impute missing values with a value that is 1 SD below the global minimum
  imputed_matrix <- assay(data_norm)
  # Replace the imputed values with the specified value
  imputed_matrix[is.na(imputed_matrix)] <- value_impute
  # Replace the assay in the original SummarizedExperiment with the imputed matrix.
  manual_imputation <- data_norm
  # Replace the assay in the original SummarizedExperiment with the imputed matrix.
  assay(manual_imputation) <- imputed_matrix


  ## Advance Imputation method ####
  ### Mixed imputation on proteins (rows)
  "(Taken from: https://bioconductor.org/packages/3.18/bioc/vignettes/DEP/inst/doc/MissingValues.html#roc-curves)
  One can also perform a mixed imputation on the proteins, which uses a MAR and MNAR imputation method on different subsets of proteins.
  First, we have to define a logical vector defining the rows that are to be imputed with the MAR method.
  Here, we consider a protein to have missing values not at random (MNAR) if it has missing values in all replicates of at least one condition."
  # Extract protein names with missing values
  # in all replicates of at least one condition

  ### Mixed Condition-Splited ####
  # 1. Create the MNAR flag table from the long-format data.
  proteins_MNAR <- get_df_long(data_norm) %>%
    mutate(across(where(is_character), as_factor)) %>%
    group_by(name, condition) %>%
    summarize(
      frac_NA = sum(is.na(intensity)) / n(), # fraction of missing
      num_NAs = sum(is.na(intensity)), # count of missing
      MNAR_flag = frac_NA >= fraction_NA, # TRUE if ≥50% missing
      .groups = "drop"
    )

  # save proteins_MNAR as xlsx
  write.xlsx(proteins_MNAR, paste0(output_path, "tables/proteins_MNAR_", name_df, ".xlsx"))


  # 2. Get the unique conditions from the SummarizedExperiment metadata.
  #    We assume that colData(data_norm)$condition exists.
  conditions <- unique(colData(data_norm)$condition)
  # Global minimum for left-censored distribution
  global_min <- min(assay(data_norm), na.rm = TRUE)
  value_impute <- global_min - sd(assay(data_norm), na.rm = TRUE) * factor_SD_impute
  # 3. Process each condition separately.
  imputed_list <- lapply(conditions, function(cond) {
    # Subset the SummarizedExperiment to only the columns (samples) for this condition.
    se_cond <- data_norm[, colData(data_norm)$condition == cond]

    # For each protein in this subset, get its MNAR flag.
    # If a protein is not found, assume it is MAR (i.e. randna = TRUE).
    # The impute function expects: TRUE = MAR (impute with knn) and FALSE = MNAR (impute as zero).
    randna_vec <- sapply(rownames(se_cond), function(prot) {
      flag <- proteins_MNAR %>%
        filter(name == prot, condition == cond) %>%
        pull(MNAR_flag)
      if (length(flag) == 0) {
        # If no record is found, assume MAR.
        return(TRUE)
      } else {
        # If MNAR_flag is TRUE (i.e. high missingness) then we want zero imputation,
        # so set randna (MAR flag) to FALSE.
        return(!flag)
      }
    })

    # Sanity check: ensure we have one flag per protein.
    if (length(randna_vec) != nrow(se_cond)) {
      stop("Mismatch in protein numbers for condition: ", cond)
    }

    # Perform the mixed imputation on this condition subset.
    imputed_se <- impute(se_cond,
      fun = "mixed",
      randna = randna_vec, # vector for this condition
      mar = "knn", # imputation method for MAR values
      mnar = "zero" # imputation method for MNAR values
    ) # imputation method for MNAR values
    return(imputed_se)
  })


  # 4. Combine the imputed data back together.
  #    Here, we assume that the protein order is the same across conditions.
  #    We combine the assays (columns) from each condition subset.
  imputed_assays <- lapply(imputed_list, assay)
  imputed_matrix <- do.call(cbind, imputed_assays)
  # value_imputation to reemplaze zeros by mannual value
  imputed_matrix[imputed_matrix == 0] <- value_impute

  # 5. Replace the assay in the original SummarizedExperiment with the imputed matrix.
  #    This creates a new SummarizedExperiment with imputed data.
  mixed_splited_imputation <- data_norm
  assay(mixed_splited_imputation) <- imputed_matrix

  proteins_MNAR <- proteins_MNAR %>%
    pivot_wider(
      names_from = condition,
      values_from = c(frac_NA, MNAR_flag, num_NAs),
      names_sep = "_"
    ) %>%
    dplyr::rename(genename = name)
  # Ahora agrega proteins_MNAR como metadata adicional:
  metadata(mixed_splited_imputation)$proteins_MNAR <- proteins_MNAR

  # Now mixed_imputation contains the imputed data.


  # Comparison SD of imputation methods ####
  compare_sd_imputations <- function(data_norm, impute_list, output_path) {
    # 1) SD before
    sd_before_df <- get_df_long(data_norm) %>%
      mutate(across(where(is.character), as_factor)) %>%
      group_by(name, condition) %>%
      summarize(
        sd_before = sd(intensity, na.rm = TRUE),
        .groups   = "drop"
      )

    # 2) SD after, for each imputed object
    sd_after_df <- bind_rows(
      lapply(names(impute_list), function(meth) {
        get_df_long(impute_list[[meth]]) %>%
          mutate(across(where(is.character), as_factor)) %>%
          group_by(name, condition) %>%
          summarize(
            sd_after = sd(intensity, na.rm = TRUE),
            .groups  = "drop"
          ) %>%
          mutate(method = meth)
      })
    )

    # 3) join before & after
    var_df <- sd_after_df %>%
      left_join(sd_before_df, by = c("name", "condition"))

    # 4) wide table of per‐protein SDs
    var_df_wide <- var_df %>%
      pivot_wider(
        id_cols = c(name, condition, sd_before),
        names_from = method,
        values_from = sd_after,
        names_prefix = "sd_after_"
      )
    assign("var_df_wide", var_df_wide, envir = .GlobalEnv)

    # export wide SD table
    dir.create(file.path(output_path, "tables"), recursive = TRUE, showWarnings = FALSE)
    write.xlsx(
      var_df_wide,
      file = file.path(output_path, "tables", "SD_before_after_imputation.xlsx"),
      overwrite = TRUE
    )
    print(var_df_wide)

    # 5) Build per‐method lm summaries and deviation metrics
    metrics_df <- var_df %>%
      filter(!is.na(sd_before), !is.na(sd_after)) %>%
      group_by(method) %>%
      nest() %>%
      mutate(
        model     = purrr::map(data, ~ lm(sd_after ~ sd_before, data = .x)),
        glance    = purrr::map(model, broom::glance),
        tidy_coef = purrr::map(model, broom::tidy)
      ) %>%
      unnest_wider(glance, names_sep = "_") %>%
      unnest(tidy_coef) %>%
      filter(term %in% c("(Intercept)", "sd_before")) %>%
      select(
        method,
        term,
        estimate,
        p.value,
        glance_r.squared,
        glance_adj.r.squared
      ) %>%
      pivot_wider(
        names_from  = term,
        values_from = c(estimate, p.value),
        names_glue  = "{term}_{.value}"
      ) %>%
      rename(
        intercept     = `(Intercept)_estimate`,
        slope         = sd_before_estimate,
        p_intercept   = `(Intercept)_p.value`,
        p_slope       = sd_before_p.value,
        r_squared     = glance_r.squared,
        adj_r_squared = glance_adj.r.squared
      ) %>%
      left_join(
        var_df %>%
          mutate(delta_sd = sd_after - sd_before) %>%
          group_by(method) %>%
          summarize(
            mean_delta = mean(delta_sd, na.rm = TRUE),
            rmse       = sqrt(mean(delta_sd^2, na.rm = TRUE)),
            .groups    = "drop"
          ),
        by = "method"
      ) %>%
      mutate(
        label = sprintf(
          "slope = %.3f (p = %.3g)\nint = %.3f (p = %.3g)\nR² = %.3f\nadjR² = %.3f\nmeanΔ = %.3f\nRMSE = %.3f",
          slope, p_slope,
          intercept, p_intercept,
          r_squared,
          adj_r_squared,
          mean_delta,
          rmse
        )
      )

    # write out the enhanced summary
    write.xlsx(
      metrics_df %>% select(-label),
      file = file.path(output_path, "tables", paste0(
        "fraction_NA_", fraction_NA,
        "_factor_SD_impute_", factor_SD_impute,
        "_SD_imputation_summary.xlsx"
      )),
      overwrite = TRUE
    )
    print(metrics_df)

    # 6) scatter + annotations
    p1 <- ggplot(var_df, aes(x = sd_before, y = sd_after)) +
      geom_point(alpha = 0.6) +
      geom_abline(
        slope = 1, intercept = 0,
        color = "red", linetype = "dashed"
      ) +
      geom_smooth(method = "lm", se = FALSE, color = "blue") +
      facet_wrap(~method, scales = "free") +
      geom_text(
        data = metrics_df,
        aes(x = Inf, y = -Inf, label = label),
        hjust = 1.1,
        vjust = -0.2,
        size = 3,
        inherit.aes = FALSE
      ) +
      labs(
        title = "SD before vs. after imputation\n(red = identity, blue = fit)",
        x     = "SD before",
        y     = "SD after"
      ) +
      theme_minimal()
    print(p1)

    # save scatter plot
    dir.create(file.path(output_path, "figures"), recursive = TRUE, showWarnings = FALSE)
    base1 <- file.path(
      output_path, "figures",
      sprintf("%02d_SD_before_after_scatter", image_number)
    )
    ggsave(paste0(base1, tiff_extension), p1, width = 8, height = 6)
    ggsave(paste0(base1, pdf_extension), p1, width = 8, height = 6)
    image_number <<- image_number + 1
  }

  # Example usage:
  imputation_list <- list(
    data_imp = data_imp,
    data_imp_man_gauss = data_imp_man_gauss,
    data_imp_knn = data_imp_knn,
    mixed_splited_imputation = mixed_splited_imputation,
    manual_imputation = manual_imputation,
    data_imp_QRILC = data_imp_QRILC
  )
  compare_sd_imputations(data_norm, imputation_list, output_path)


  ## Test for differential analysis for imputation methods
  # We perform differential analysis on the different imputated data sets. The following datasets are compared: No imputation, knn imputation, MinProb imputation, Mixed imputation.

  # Exporting results for data analysis ####

  save(data_imp, file = paste0(output_path, "RData/data_imp_", name_df, ".RData"))
  save(data_imp_man_gauss, file = paste0(output_path, "RData/data_imp_man_gauss_", name_df, ".RData"))
  save(data_imp_knn, file = paste0(output_path, "RData/data_imp_knn_", name_df, ".RData"))
  save(mixed_splited_imputation, file = paste0(output_path, "RData/mixed_splited_imputation_", name_df, ".RData"))
  save(manual_imputation, file = paste0(output_path, "RData/manual_imputation_", name_df, ".RData"))
  save(data_imp_QRILC, file = paste0(output_path, "RData/data_imp_QRILC_", name_df, ".RData"))

  load(paste0(output_path, "RData/data_imp_", name_df, ".RData"))
  load(paste0(output_path, "RData/data_imp_man_gauss_", name_df, ".RData"))
  load(paste0(output_path, "RData/data_imp_knn_", name_df, ".RData"))
  load(paste0(output_path, "RData/manual_imputation_", name_df, ".RData"))
  load(paste0(output_path, "RData/mixed_splited_imputation_", name_df, ".RData"))
  load(paste0(output_path, "RData/data_norm_", name_df, ".RData"))
  load(paste0(output_path, "RData/data_imp_QRILC_", name_df, ".RData"))

  # Dynamically assign the data frame to the global environment
  assign("data_imp", data_imp, envir = .GlobalEnv)
  assign("data_imp_man_gauss", data_imp_man_gauss, envir = .GlobalEnv)
  assign("data_imp_knn", data_imp_knn, envir = .GlobalEnv)
  assign("manual_imputation", manual_imputation, envir = .GlobalEnv)
  assign("mixed_splited_imputation", mixed_splited_imputation, envir = .GlobalEnv)
  assign("data_imp_QRILC", data_imp_QRILC, envir = .GlobalEnv)

  # Visualization of imputation effect
  p <- plot_imputation(data_norm, data_imp, data_imp_man_gauss, data_imp_knn, manual_imputation, mixed_splited_imputation, data_imp_QRILC)
  print(p)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_protein_imputation_distribution_", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), p, width = 8, height = 6, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), p, width = 8, height = 6, units = "in")
  image_number <<- image_number + 1
  # For our dataset, knn and mixed imputation result in less identified differential expressed proteins compared to the no imputation and MinProb. No imputation results in the identification of the most differentially expressed proteins in our dataset with many proteins missing values.
  # Note that the performance of the different imputation methods is data set-dependent. It is recommended to always carefully check the effect of filtering and data imputation on your results.

  # PCA plots ####
  dep_analysis_norm <- analyze_dep(data_norm, type = "manual", control = NULL, alpha = p_val, lfc = FC, test = comparisons, design_formula = formula(~ 0 + condition))
  dep_analysis_min <- analyze_dep(data_imp, type = "manual", control = NULL, alpha = p_val, lfc = FC, test = comparisons, design_formula = formula(~ 0 + condition))
  dep_analysis_manual_gauss <- analyze_dep(data_imp_man_gauss, type = "manual", control = NULL, alpha = p_val, lfc = FC, test = comparisons, design_formula = formula(~ 0 + condition))
  dep_analysis_knn <- analyze_dep(data_imp_knn, type = "manual", control = NULL, alpha = p_val, lfc = FC, test = comparisons, design_formula = formula(~ 0 + condition))
  dep_analysis_manual_value <- analyze_dep(manual_imputation, type = "manual", control = NULL, alpha = p_val, lfc = FC, test = comparisons, design_formula = formula(~ 0 + condition))
  dep_analysis_mixed_splited <- analyze_dep(mixed_splited_imputation, type = "manual", control = NULL, alpha = p_val, lfc = FC, test = comparisons, design_formula = formula(~ 0 + condition))
  dep_analysis_QRILC <- analyze_dep(data_imp_QRILC, type = "manual", control = NULL, alpha = p_val, lfc = FC, test = comparisons, design_formula = formula(~ 0 + condition))


  # PCA MinProb
  pca_minprob <- plot_pca(dep_analysis_min, x = 2, y = 1, n = (length(dep_analysis_min)), point_size = 4) + ggtitle("PCA MinProb", subtitle = paste0(length(dep_analysis_min), " variable proteins"))
  print(pca_minprob)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_PCA_MinProb_", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), pca_minprob, width = 8, height = 6, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), pca_minprob, width = 8, height = 6, units = "in")
  image_number <<- image_number + 1

  # PCA Manual
  pca_manual_gauss <- plot_pca(dep_analysis_manual_gauss, x = 2, y = 1, n = (length(dep_analysis_manual_gauss)), point_size = 4) + ggtitle("PCA Manual", subtitle = paste0(length(dep_analysis_manual_gauss), " variable proteins"))
  print(pca_manual_gauss)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_PCA_Manual_Gauss", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), pca_manual_gauss, width = 8, height = 6, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), pca_manual_gauss, width = 8, height = 6, units = "in")
  image_number <<- image_number + 1

  # PCA KNN
  pca_knn <- plot_pca(dep_analysis_knn, x = 2, y = 1, n = (length(dep_analysis_knn)), point_size = 4) + ggtitle("PCA KNN", subtitle = paste0(length(dep_analysis_knn), " variable proteins"))
  print(pca_knn)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_PCA_KNN_", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), pca_knn, width = 8, height = 6, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), pca_knn, width = 8, height = 6, units = "in")
  image_number <<- image_number + 1

  # PCA Mixed Splited
  pca_mixed_splited <- plot_pca(dep_analysis_mixed_splited, x = 2, y = 1, n = (length(dep_analysis_mixed_splited)), point_size = 4) + ggtitle("PCA Mixed Splited Imputation", subtitle = paste0(length(dep_analysis_mixed_splited), " variable proteins"))
  print(pca_mixed_splited)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_PCA_Splited_Mixed_", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), pca_mixed_splited, width = 8, height = 6, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), pca_mixed_splited, width = 8, height = 6, units = "in")
  image_number <<- image_number + 1

  # PCA manual_imputation
  pca_manual_value <- plot_pca(dep_analysis_manual_value, x = 2, y = 1, n = (length(dep_analysis_manual_value)), point_size = 4) + ggtitle("PCA Manual Imputation", subtitle = paste0(length(dep_analysis_manual_value), " variable proteins"))
  print(pca_manual_value)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_PCA_manual_imputation_", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), pca_manual_value, width = 8, height = 6, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), pca_manual_value, width = 8, height = 6, units = "in")
  image_number <<- image_number + 1

  # PCA QRILC
  pca_QRILC <- plot_pca(dep_analysis_QRILC, x = 2, y = 1, n = (length(dep_analysis_QRILC)), point_size = 4) + ggtitle("PCA QRILC", subtitle = paste0(length(dep_analysis_QRILC), " variable proteins"))
  print(pca_QRILC)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_PCA_QRILC_", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), pca_QRILC, width = 8, height = 6, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), pca_QRILC, width = 8, height = 6, units = "in")
  image_number <<- image_number + 1
  # Patchwork of PCA plots
  pca_combined <- (pca_minprob | pca_manual_gauss) / (pca_manual_value | pca_QRILC) / (pca_knn | pca_mixed_splited) +
    plot_annotation(
      title = "PCA plots of different imputation methods",
      # Label each plot
      theme = theme(plot.title = element_text(size = 16, hjust = 0.5)),
      tag_levels = "A"
    )
  print(pca_combined)
  filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_PCA_imputation_methods_", name_df) # Name of output file.
  ggsave(paste0(filename, tiff_extension), pca_combined, width = 12, height = 18, units = "in") # Adjust size according to your needs.
  ggsave(paste0(filename, pdf_extension), pca_combined, width = 12, height = 18, units = "in")
  image_number <<- image_number + 1


  # PCA-plot Outlier identification ####
  # Function to identify outliers within a group using robust PCA with error handling
  run_pca_hubert_analysis <- function(pca_list, k = 5) {
    all_outliers <- data.frame(Sample = pca_list[[1]]$data$rowname)
    results <- list()
    plots <- list() # Store plots for patchwork

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

      # Generar gráfico y almacenarlo
      plot_df <- distances_df
      p <- ggplot(plot_df, aes(x = Score_Distance, y = Orthogonal_Distance, label = Sample)) +
        geom_point(aes(color = Outlier), size = 3) +
        geom_text(hjust = 1.2, vjust = 0.5, size = 3) +
        scale_color_manual(values = c("black", "red")) +
        theme_minimal() +
        labs(
          title = paste("DD-plot PCAHubert:", pca_name),
          subtitle = "Outliers in red",
          x = "Score Distance",
          y = "Orthogonal Distance"
        )

      # Store plot in list instead of printing
      plots[[pca_name]] <- p
      results[[pca_name]] <- distances_df
    }

    # Create combined plot using patchwork
    if (length(plots) > 0) {
      # Combine all plots using patchwork
      combined_plot <- wrap_plots(plots, ncol = 2) # Adjust ncol as needed

      # Save combined plot
      filename <- paste0(output_path, "figures/", sprintf("%02d", image_number), "_Combined_DD_Plots")
      ggsave(paste0(filename, ".tiff"), combined_plot, width = 16, height = 12, units = "in")
      ggsave(paste0(filename, ".pdf"), combined_plot, width = 16, height = 12, units = "in")

      # Print combined plot
      print(combined_plot)

      image_number <<- image_number + 1
    }

    # Guardar tabla consolidada de outliers
    write.xlsx(all_outliers, paste0(output_path, "tables/PCA_based_outliers.xlsx"))
    # Save results
    write.xlsx(results, paste0(output_path, "tables/PCA_detailed_outlier_results.xlsx"))
    return(list("outlier_summary" = all_outliers, "detailed_results" = results, "combined_plot" = combined_plot))
  }

  # List of PCA objects
  pca_objects <- list(
    pca_minprob = pca_minprob,
    pca_manual_gauss = pca_manual_gauss,
    pca_knn = pca_knn,
    pca_manual_value = pca_manual_value,
    pca_mixed_splited = pca_mixed_splited,
    pca_QRILC = pca_QRILC
  )

  # Run analysis
  pca_results <- run_pca_hubert_analysis(pca_objects)
  print(pca_results$outlier_summary)


  save(pca_results, file = paste0(output_path, "RData/pca_outliers_results.RData"))


  # Differential Expression ####
  # Function that wraps around test_diff, add_rejections and get_results functions
  DE_analysis <- function(se) {
    se %>%
      test_diff(., type = "manual", test = comparisons) %>%
      add_rejections(., alpha = p_val, lfc = log2(FC)) %>%
      get_results()
  }

  # DE analysis on no, knn, MinProb and mixed imputation
  no_imputation_results <- DE_analysis(data_norm)
  MinProb_imputation_results <- DE_analysis(data_imp)
  Manual_imputation_results <- DE_analysis(data_imp_man_gauss)
  knn_imputation_results <- DE_analysis(data_imp_knn)
  manual_imputation_results <- DE_analysis(manual_imputation)
  mixed_splited_imputation_results <- DE_analysis(mixed_splited_imputation)
  QRILC_imputation_results <- DE_analysis(data_imp_QRILC)


  ### Number of identified differentially expressed proteins
  # As an initial parameter we look at the number of differentially expressed proteins identified (adjusted P ≤ 0.05 and Fold-Change > 1).

  # Function to extract number of DE proteins
  DE_prots <- function(results) {
    prots <- tibble(
      Dataset = gsub("_results", "", results),
      significant_proteins = get(results) %>%
        filter(significant)
    )
    number_prots <- tibble(
      Dataset = gsub("_results", "", results),
      significant_proteins = get(results) %>%
        filter(significant) %>%
        nrow()
    )
  }

  # Number of significant proteins
  objects <- c(
    "no_imputation_results",
    "MinProb_imputation_results",
    "Manual_imputation_results",
    "knn_imputation_results",
    "manual_imputation_results",
    "mixed_splited_imputation_results",
    "QRILC_imputation_results"
  )

  dep_imputed_prots <- map_df(objects, DE_prots)
  # dep_imputed_prots$imputation_mnar <- mnar_var
  dep_imputed_prots$global_min <- global_min
  dep_imputed_prots$imputation_value <- value_impute

  print(dep_imputed_prots)
  dep_imputed_prots %>%
    write_xlsx(paste0(output_path, "tables/dep_imputed_prots_", name_df, ".xlsx"))
}
