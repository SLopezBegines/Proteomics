#!/usr/bin/env Rscript
# Load the protein table, filter, normalise (vsn) and impute missing values.
# Wraps code/01_loading_data.R and code/03_cleaning_data_mixed_imputation.R.
this_dir <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])))
source(file.path(this_dir, "pipeline_common.R"))

a <- parse_args()
req <- c("protein_groups", "design", "contrasts", "dataset", "organism", "state", "fraction_na",
         "factor_sd_impute", "factor_sd_scale", "mnar_var", "p_val", "p_val_low", "fc")
for (k in req) if (is.null(a[[k]])) die("Missing --", k)
a$state <- normalizePath(a$state, mustWork = TRUE)
pg <- normalizePath(a$protein_groups, mustWork = TRUE)
design <- read_design(a$design)
ctx <- list(
  dataset = a$dataset, organism = a$organism,
  exp_design = build_exp_design(design, a$dataset), comparisons = read_contrasts(a$contrasts),
  p_val = as.numeric(a$p_val), p_val_low = as.numeric(a$p_val_low), fc = as.numeric(a$fc),
  image_number = 1
)
args_num <- list(fraction_NA = as.numeric(a$fraction_na), factor_SD_impute = as.numeric(a$factor_sd_impute),
                 factor_SD_scale = as.numeric(a$factor_sd_scale))

prot_data <- read_protein_table(pg)          # read BEFORE setup_pipeline() changes the working directory
setup_pipeline(ctx, a$state)
Exp_design <- ctx$exp_design
comparisons <- ctx$comparisons

source("../code/01_loading_data.R")
n_in <- nrow(prot_data)
loading_data(prot_data)                       # writes tables/*_selected_data|contaminants|cleaned_data.xlsx
cat(sprintf("Proteins: %d in input, %d after contaminant/NA-gene filtering\n", n_in, nrow(prot_data)))

source("../code/03_cleaning_data_mixed_imputation.R")
image_number <- 1
data_cleaning(prot_data, Exp_design, comparisons, output_path,
              fraction_NA = args_num$fraction_NA, factor_SD_scale = args_num$factor_SD_scale,
              factor_SD_impute = args_num$factor_SD_impute, mnar_var = a$mnar_var)

save_context(ctx, a$state)
write_versions("PREPROCESS_QC")
