# Global Variables ####

# Setting threshold for p-value and Fold-Change
p_val <- 0.05
p_val_low <- 0.01
FC <- 0.5

# Global variables
tiff_extension <- ".tiff"
pdf_extension <- ".pdf" # Vectorial format
# kegg_organism = "dre"#dre for Danio rerio, hsa for Homo sapiens, mmu for Mus musculus
# species <- 7955 #9606 for Human, 10090 for mouse, 7955 for zebrafish
# organism <- "org.Dr.eg.db" #"org.Hs.eg.db"
keyType <- "UNIPROT"
KEGGkeyType <- "uniprot"

# Color by condition
# Generate 4 colors from the Set2 palette
# Set number of color on basis of number of conditions
my_colors <- setNames(RColorBrewer::brewer.pal(4, "Set2"), c("CTRL", "CTRL_K", "MDD", "MDD_K"))

# Creat directories to store data
create_directories <- function(base_path) {
  # Helper function to create a directory if it does not exist
  create_dir_if_not_exists <- function(path) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
      cat("Directory created:", path, "\n")
    } else {
      # cat("Directory already exists:", path, "\n")
    }
  }

  # Create the base directory and its subdirectories
  create_dir_if_not_exists(base_path)
  create_dir_if_not_exists(file.path(base_path, "tables"))
  create_dir_if_not_exists(file.path(base_path, "figures"))
  create_dir_if_not_exists(file.path(base_path, "RData"))
  create_dir_if_not_exists(file.path(base_path, "figures", "VennDiagram"))
  create_dir_if_not_exists(file.path(base_path, "figures", "summary"))
  create_dir_if_not_exists(file.path(base_path, "figures", "enrichGO"))
  create_dir_if_not_exists(file.path(base_path, "figures", "gseGO"))
  create_dir_if_not_exists(file.path(base_path, "figures", "KEGG_GO"))
  create_dir_if_not_exists(file.path(base_path, "figures", "GO_adj"))
  create_dir_if_not_exists(file.path(base_path, "figures", "GO_adj"))
  create_dir_if_not_exists(file.path(base_path, "figures", "gseGO_adj"))
  create_dir_if_not_exists(file.path(base_path, "figures", "gseKEGG_adj"))
  create_dir_if_not_exists(file.path(base_path, "figures", "panther"))
  create_dir_if_not_exists(file.path(base_path, "figures", "rbioapi"))
}
create_directories(output_path)


# add prefix lfq_ to LFQ imputed columns
