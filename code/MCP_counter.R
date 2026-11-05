# Loading libraries ####

# install.packages("devtools")
# library(devtools)
# devtools::install_github("cit-bioinfo/mMCP-counter")
# renv::snapshot(lockfile = "renv-R543.lockº)


library("mMCPcounter")
# load data ####

load("/mnt/datos/github_repos/propietary/Proteomics-dep/mains/output/joan_453/RData/dep_analysis.RData")
data_results <- readxl::read_xlsx(path = file.path("mains", output_path, "tables", "data_results.xlsx"))
data <- data_results |>
  dplyr::select(starts_with("lfq"))
data_matrix <- as.matrix(data)
row.names(data_matrix) <- data_results$name

# MCP-counter estimation ####
mcp_data <- mMCPcounter.estimate(data_matrix, features = "Gene.Symbol", genomeVersion = "GCRm38")

# Convert to data frame and add ID column
mcp_data <- as.data.frame(mcp_data)
mcp_data <- mcp_data |> mutate(ID = row.names(mcp_data))

#  detele rownames
mcp_data <- mcp_data |>
  remove_rownames() |>
  dplyr::select(ID, everything())

# Convert to long format
mcp_data_long <- mcp_data |>
  pivot_longer(cols = starts_with("lfq"), names_to = "cols", values_to = "value") |>
  mutate(
    group = stringr::str_extract(cols, "(WT|HT|HM)_(s|c)"),
    replica = as.factor(stringr::str_extract(cols, "\\d+")),
    tissue = case_when(
      stringr::str_detect(group, "_s$") ~ "skin",
      stringr::str_detect(group, "_c$") ~ "colon",
      .default = NA_character_
    )
  )

# reorder levels
mcp_data_long$tissue <- factor(mcp_data_long$tissue, c("skin", "colon"))
mcp_data_long$group <- factor(mcp_data_long$group, c("WT_s", "HT_s", "HM_s", "WT_c", "HT_c", "HM_c"))

# MCP Means ####
mcp_mean <- mcp_data_long |>
  dplyr::group_by(ID, group) |>
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE),
    n = n(),
    sem = sd_value / sqrt(n),
    .groups = "drop"
  )

# Statistical analys ####
# Compare groups using ANOVA and Tukey's HSD test for each ID, tissue and group

## split by tissue ####
# skin and colon are separate experiments: one list element per tissue,
# dropping unused group levels so each ANOVA only sees its own 3 groups
mcp_by_tissue <- mcp_data_long |>
  dplyr::filter(!is.na(value)) |>
  (\(df) split(df, df$tissue))() |>
  purrr::map(\(df) dplyr::mutate(df, group = droplevels(group)))

## assumptions ####
# Shapiro on ANOVA residuals; Levene on value ~ group (centres each group, i.e. works on the model residuals)
mcp_assumptions <- purrr::imap_dfr(mcp_by_tissue, \(df, tis) {
  shapiro <- df |>
    dplyr::group_by(ID) |>
    dplyr::mutate(resid = resid(lm(value ~ group))) |>
    rstatix::shapiro_test(resid) |>
    dplyr::select(ID, shapiro_p = p)
  levene <- df |>
    dplyr::group_by(ID) |>
    rstatix::levene_test(value ~ group) |>
    dplyr::select(ID, levene_p = p)
  dplyr::left_join(shapiro, levene, by = "ID") |>
    dplyr::mutate(tissue = tis, .before = 1)
})

## ANOVA: value ~ group, per ID within each tissue ####
mcp_anova <- purrr::imap_dfr(mcp_by_tissue, \(df, tis) {
  df |>
    dplyr::group_by(ID) |>
    rstatix::anova_test(value ~ group) |>
    as.data.frame() |> # drop grouped_anova_test class, incompatible with bind_rows
    dplyr::mutate(tissue = tis, .before = 1)
}) |>
  dplyr::group_by(tissue) |>
  dplyr::mutate(p.adj.fdr = p.adjust(p, method = "BH")) |> # FDR across IDs within a tissue
  dplyr::ungroup() |>
  rstatix::add_significance(p.col = "p.adj.fdr")

## Tukey HSD post-hoc ####
mcp_tukey <- purrr::imap_dfr(mcp_by_tissue, \(df, tis) {
  df |>
    dplyr::group_by(ID) |>
    rstatix::tukey_hsd(value ~ group) |>
    as.data.frame() |>
    dplyr::mutate(tissue = tis, .before = 1)
})

### save tables ####
writexl::write_xlsx(
  list(anova = mcp_anova, tukey = mcp_tukey, assumptions = mcp_assumptions),
  path = file.path("mains", output_path, "tables", "mcp_counter_stats.xlsx")
)


# Visualization ####
#

## plots with significance ####
# Tukey brackets placed above each ID's data; only significant comparisons are drawn
add_brackets <- function(tissue_name, plot) {
  ypos <- mcp_by_tissue[[tissue_name]] |>
    dplyr::group_by(ID) |>
    dplyr::summarise(ymax = max(value), yrange = diff(range(value)), .groups = "drop")
  brackets <- mcp_tukey |>
    dplyr::filter(tissue == tissue_name, p.adj < 0.05) |>
    dplyr::left_join(ypos, by = "ID") |>
    dplyr::group_by(ID) |>
    dplyr::mutate(y.position = ymax + yrange * 0.12 * dplyr::row_number()) |>
    dplyr::ungroup()
  plot +
    facet_wrap(vars(ID), scales = "free_y") +
    ggpubr::stat_pvalue_manual(brackets, label = "p.adj.signif", tip.length = 0.01, inherit.aes = FALSE)
}
## skin data ####
p_skin <- mcp_data_long |>
  dplyr::filter(tissue == "skin") |>
  ggplot() +
  aes(x = group, y = value, fill = ID) +
  geom_jitter(aes(color = replica)) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(vars(ID))

## colon data ####
p_colon <- mcp_data_long |>
  dplyr::filter(tissue == "colon") |>
  ggplot() +
  aes(x = group, y = value, fill = ID) +
  geom_jitter(aes(color = replica)) +
  geom_boxplot() +
  scale_fill_hue(direction = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(vars(ID))

p_skin_sig <- add_brackets("skin", p_skin)
p_colon_sig <- add_brackets("colon", p_colon)

### save plots ####
save_plot("plot_MCP_skin_sig", p_skin_sig, subdir = "mcp_counter")
save_plot("plot_MCP_colon_sig", p_colon_sig, subdir = "mcp_counter")


# heatmap ####
# heatmap of mean values, with dendrograms for ID and group
heatmap_data <- mcp_mean |>
  dplyr::select(ID, group, mean_value) |>
  tidyr::pivot_wider(names_from = group, values_from = mean_value) |>
  tibble::column_to_rownames("ID") |>
  as.matrix()

pheatmap <- pheatmap::pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  color = viridis::viridis(100),
  fontsize_row = 8,
  fontsize_col = 10,
  main = "MCP-counter mean values"
)
