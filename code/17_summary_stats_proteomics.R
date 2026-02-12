# Proteomic stats


# Load data results
data_brain_CLN3 <- readxl::read_xlsx("mains/output/CLN3/Brains/tables/data_results.xlsx")
data_larvae_CLN3 <- readxl::read_xlsx("mains/output/CLN3/Larvae/outliers_removed/KNN/tables/data_results.xlsx")

data_brain_CLN12 <- readxl::read_xlsx("mains/output/CLN12/Brains/outlier_removed/tables/data_results.xlsx")
data_larvae_CLN12 <- readxl::read_xlsx("mains/output/CLN12/Larvae/outlier_removed/tables/data_results.xlsx")
# Define your comparisons
comparisons <- c("CLN3_Lux1_vs_WT", "CLN3_Lux2_vs_WT")
# If you want to include CLN12 comparisons, uncomment the following line:
comparisons <- c("CLN12_vs_WT")
# Create a named list of your data frames
list_results <- list(data_brain_CLN3, data_larvae_CLN3)
names(list_results) <- c("Brain_CLN3", "Larvae_CLN3")
# If you want to include CLN12 data, uncomment the following lines:
list_results <- list(data_brain_CLN12, data_larvae_CLN12)
names(list_results) <- c("Brain_CLN12", "Larvae_CLN12")

# Number of identified proteins ####
# Function to calculate proteomic statistics
calculate_proteomic_stats <- function(data, comparisons, dataset_name = NULL) {
  # Initialize results list
  stats_list <- list()

  # Add total proteins
  stats_list[[1]] <- tibble(
    Dataset = ifelse(is.null(dataset_name), "Unknown", dataset_name),
    Metric = "Total proteins",
    Comparison = "All",
    Count = nrow(data),
    Percentage = 100
  )

  # Loop through each comparison
  for (comp in comparisons) {
    # Column names
    pval_col <- paste0(comp, "_p.val")
    padj_col <- paste0(comp, "_p.adj")
    ratio_col <- paste0(comp, "_ratio")

    # Check if columns exist
    if (!all(c(pval_col, padj_col, ratio_col) %in% colnames(data))) {
      warning(paste("Columns for comparison", comp, "not found in dataset", dataset_name))
      next
    }

    # Calculate statistics for this comparison
    comp_stats <- tibble(
      Dataset = ifelse(is.null(dataset_name), "Unknown", dataset_name),
      Metric = c(
        "p.val < 0.05",
        "p.adj < 0.05",
        "p.val < 0.05 & Ratio > 0 (upregulated)",
        "p.val < 0.05 & Ratio < 0 (downregulated)",
        "p.adj < 0.05 & Ratio > 0 (upregulated)",
        "p.adj < 0.05 & Ratio < 0 (downregulated)"
      ),
      Comparison = comp,
      Count = c(
        sum(data[[pval_col]] < 0.05, na.rm = TRUE),
        sum(data[[padj_col]] < 0.05, na.rm = TRUE),
        sum(data[[pval_col]] < 0.05 & data[[ratio_col]] > 0, na.rm = TRUE),
        sum(data[[pval_col]] < 0.05 & data[[ratio_col]] < 0, na.rm = TRUE),
        sum(data[[padj_col]] < 0.05 & data[[ratio_col]] > 0, na.rm = TRUE),
        sum(data[[padj_col]] < 0.05 & data[[ratio_col]] < 0, na.rm = TRUE)
      )
    ) %>%
      mutate(Percentage = round((Count / nrow(data)) * 100, 2))

    stats_list[[length(stats_list) + 1]] <- comp_stats
  }

  # Combine all statistics
  stats_table <- bind_rows(stats_list)

  return(stats_table)
}

# Apply function to a list and maintain names (combined output)
apply_stats_to_list <- function(data_list, comparisons) {
  # Apply function to each element with its name
  stats_results <- imap_dfr(data_list, ~ calculate_proteomic_stats(.x, comparisons, .y))

  return(stats_results)
}

# Apply function to a list and maintain names (separate outputs)
apply_stats_to_list_separate <- function(data_list, comparisons) {
  # Apply function to each element with its name
  stats_results <- imap(data_list, ~ calculate_proteomic_stats(.x, comparisons, .y))

  return(stats_results)
}

# Example usage:
# Define your comparisons
# comparisons <- c("CLN3_Lux1_vs_WT", "CLN3_Lux2_vs_WT")
#
# Create a named list of your data frames
# data_list <- list(
#   experiment1 = read_excel("data1.xlsx"),
#   experiment2 = read_excel("data2.xlsx")
# )
#
# Option 1: Combined table with all datasets
all_stats <- apply_stats_to_list(list_results, comparisons)
print(all_stats)
# save as xlsx
writexl::write_xlsx(all_stats, "CLN3_all_proteomic_stats.xlsx")
# write_csv(all_stats, "all_proteomic_stats.csv")
#
# Option 2: Separate tables per dataset (maintains list structure)
stats_list <- apply_stats_to_list_separate(list_results, comparisons)
print(stats_list$experiment1)
#
# Example with different number of comparisons:
# comparisons_extended <- c("CLN3_Lux1_vs_WT", "CLN3_Lux2_vs_WT", "CLN3_Lux3_vs_WT")
# all_stats_extended <- apply_stats_to_list(data_list, comparisons_extended)

# Plots ####

# Function to plot statistics results
plot_proteomic_stats <- function(stats_table, plot_type = "all") {
  # Filter out total proteins for most plots
  stats_filtered <- stats_table %>%
    filter(Metric != "Total proteins")

  if (plot_type == "all" || plot_type == "summary") {
    # Plot 1: Summary of significant proteins (p.val and p.adj < 0.05)
    p1 <- stats_filtered %>%
      filter(Metric %in% c("p.val < 0.05", "p.adj < 0.05")) %>%
      ggplot(aes(x = Comparison, y = Count, fill = Metric)) +
      geom_col(position = "dodge", alpha = 0.8) +
      geom_text(aes(label = Count),
        position = position_dodge(width = 0.9),
        vjust = -0.5, size = 3
      ) +
      facet_wrap(~Dataset, scales = "free_x") +
      labs(
        title = "Significant Proteins by Comparison",
        subtitle = "Number of proteins with p-value or adjusted p-value < 0.05",
        x = "Comparison", y = "Number of Proteins", fill = "Threshold"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )

    print(p1)
  }

  if (plot_type == "all" || plot_type == "regulation") {
    # Plot 2: Up/Down regulation (p.val < 0.05)
    p2 <- stats_filtered %>%
      filter(grepl("p.val.*Ratio", Metric)) %>%
      mutate(
        Direction = ifelse(grepl("upregulated", Metric), "Upregulated", "Downregulated"),
        Direction = factor(Direction, levels = c("Upregulated", "Downregulated"))
      ) %>%
      ggplot(aes(x = Comparison, y = Count, fill = Direction)) +
      geom_col(position = "dodge", alpha = 0.8) +
      geom_text(aes(label = Count),
        position = position_dodge(width = 0.9),
        vjust = -0.5, size = 3
      ) +
      facet_wrap(~Dataset, scales = "free_x") +
      scale_fill_manual(values = c("Upregulated" = "#E74C3C", "Downregulated" = "#3498DB")) +
      labs(
        title = "Protein Regulation Direction (p-value < 0.05)",
        subtitle = "Number of upregulated vs downregulated proteins",
        x = "Comparison", y = "Number of Proteins", fill = "Regulation"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )

    print(p2)
  }

  if (plot_type == "all" || plot_type == "regulation_adj") {
    # Plot 3: Up/Down regulation (p.adj < 0.05)
    p3 <- stats_filtered %>%
      filter(grepl("p.adj.*Ratio", Metric)) %>%
      mutate(
        Direction = ifelse(grepl("upregulated", Metric), "Upregulated", "Downregulated"),
        Direction = factor(Direction, levels = c("Upregulated", "Downregulated"))
      ) %>%
      ggplot(aes(x = Comparison, y = Count, fill = Direction)) +
      geom_col(position = "dodge", alpha = 0.8) +
      geom_text(aes(label = Count),
        position = position_dodge(width = 0.9),
        vjust = -0.5, size = 3
      ) +
      facet_wrap(~Dataset, scales = "free_x") +
      scale_fill_manual(values = c("Upregulated" = "#E74C3C", "Downregulated" = "#3498DB")) +
      labs(
        title = "Protein Regulation Direction (adjusted p-value < 0.05)",
        subtitle = "Number of upregulated vs downregulated proteins (FDR corrected)",
        x = "Comparison", y = "Number of Proteins", fill = "Regulation"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )

    print(p3)
  }

  if (plot_type == "all" || plot_type == "heatmap") {
    # Plot 4: Heatmap of counts
    p4 <- stats_filtered %>%
      mutate(Metric_short = case_when(
        Metric == "p.val < 0.05" ~ "p-val sig",
        Metric == "p.adj < 0.05" ~ "p-adj sig",
        grepl("p.val.*upregulated", Metric) ~ "p-val UP",
        grepl("p.val.*downregulated", Metric) ~ "p-val DOWN",
        grepl("p.adj.*upregulated", Metric) ~ "p-adj UP",
        grepl("p.adj.*downregulated", Metric) ~ "p-adj DOWN",
        TRUE ~ Metric
      )) %>%
      ggplot(aes(x = Comparison, y = Metric_short, fill = Count)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = Count), color = "white", fontface = "bold") +
      facet_wrap(~Dataset, scales = "free_x") +
      scale_fill_gradient(low = "#3498DB", high = "#E74C3C") +
      labs(
        title = "Proteomic Statistics Heatmap",
        subtitle = "Overview of all metrics across comparisons",
        x = "Comparison", y = "Metric", fill = "Count"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )

    print(p4)
  }

  if (plot_type == "all" || plot_type == "comparison") {
    # Plot 5: Comparison between datasets (if multiple datasets)
    if (length(unique(stats_table$Dataset)) > 1) {
      p5 <- stats_filtered %>%
        filter(Metric %in% c("p.adj < 0.05")) %>%
        ggplot(aes(x = Dataset, y = Count, fill = Comparison)) +
        geom_col(position = "dodge", alpha = 0.8) +
        geom_text(aes(label = Count),
          position = position_dodge(width = 0.9),
          vjust = -0.5, size = 3
        ) +
        labs(
          title = "Dataset Comparison",
          subtitle = "Significant proteins (adjusted p-value < 0.05) across datasets",
          x = "Dataset", y = "Number of Proteins", fill = "Comparison"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom"
        )

      print(p5)
    }
  }
}



# Example usage:
# all_stats <- apply_stats_to_list(data_list, comparisons)
#
# # Plot all visualizations
plot_proteomic_stats(all_stats, plot_type = "all")
#
# # Or plot specific types:
plot_proteomic_stats(all_stats, plot_type = "summary")
plot_proteomic_stats(all_stats, plot_type = "regulation")
plot_proteomic_stats(all_stats, plot_type = "regulation_adj")
plot_proteomic_stats(all_stats, plot_type = "heatmap")
plot_proteomic_stats(all_stats, plot_type = "comparison")
