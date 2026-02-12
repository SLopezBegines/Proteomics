

# libraries ####

source("../code/00_packages.R")


# Load data ####

#Load dep_analysis object
load(paste0(output_path,"RData/dep_analysis.RData"))
#Load data results tables
data_results <- readxl::read_xlsx(path = paste0(output_path,"tables/data_results.xlsx"))
significative_data <- readxl::read_xlsx(path = paste0(output_path,"tables/significative_data.xlsx"))
sig_adjusted_data <- readxl::read_xlsx(path = paste0(output_path,"tables/sig_adjusted_data.xlsx"))


# Plots ####
## PCA plot ####
p <- plot_pca(dep_analysis,x = 2, y = 1, n = (length(dep_analysis)), point_size = 4)
print(p)
filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_PCA")# Name of output file. 
ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Adjust size according to your needs.
ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  
image_number <<- image_number + 1

## PCA plot loadings ####
# get LFQ column numbers
LFQ_columns <- data_results %>% 
  dplyr::select(starts_with("lfq_")) %>% 
  colnames()


pca_matrix <- data_results[,c("significant", "significance",LFQ_columns)]
pca_matrix <- as.matrix(pca_matrix)
row.names(pca_matrix) <- data_results$name

pca_res <- stats::prcomp(pca_matrix, scale. = FALSE)
#pca_res$group       <- c(rep("WT", 4), rep("CLN3_Lux1", 4), rep("CLN3_Lux2", 4), rep("CLN12",4))
pca_res$repetition <- rep(c(rep("1", 1), rep("2", 1),rep("3", 1),rep("4", 1)), 4)

#loadings by adj_p_value
pca <- autoplot(pca_res, data = data_results, colour= "significant",shape = TRUE, 
                label = FALSE, label.size = 3,
                loadings = FALSE, loadings.colour = "blue",
                loadings.label = TRUE, loadings.label.size = 3, frame = FALSE, frame.type = "norm")
print(pca)
filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_PCA_loadings_p_adj")  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
ggsave(paste0(filename,pdf_extension), pca, width = 5, height = 5, units = "in")  # Vectorial format
ggsave(paste0(filename,tiff_extension), pca, width = 5, height = 5, units = "in")  # Tiff format
image_number <<- image_number + 1
#loadings by p_value
pca <- autoplot(pca_res, data = data_results, colour= "significance",shape = TRUE, 
                label = FALSE, label.size = 3,
                loadings = FALSE, loadings.colour = "blue",
                loadings.label = TRUE, loadings.label.size = 3, frame = FALSE, frame.type = "norm")
print(pca)
filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_PCA_loadings_p_value")  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
ggsave(paste0(filename,pdf_extension), pca, width = 5, height = 5, units = "in")  # Vectorial format
ggsave(paste0(filename,tiff_extension), pca, width = 5, height = 5, units = "in")  # Tiff format
image_number <<- image_number + 1



## Heatmaps ####
### Correlation matrix ####

# Plot the Pearson correlation matrix
# Step 1: Call the pdf command to start the plot
filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_Correlation_matrix")# Name of output file.
pdf(file = paste0(filename,pdf_extension)) # The height of the plot in inches
p <- DEP::plot_cor(dep_analysis, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
print(p)
dev.off()
tiff(file = paste0(filename, tiff_extension)) # The height of the plot in inches
print(p)
dev.off()
image_number <<- image_number + 1
### HeatMap ####
filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_Heatmap_significant")# Name of output file.
pdf(file = paste0(filename,pdf_extension)) 
p <- DEP::plot_heatmap(dep_analysis, type = "centered",  kmeans = TRUE, 
                  k = 2, col_limit = 2, show_row_names = TRUE,
                  indicate = c("condition", "replicate"))
print(p)
dev.off()
tiff(file = paste0(filename,tiff_extension)) # The height of the plot in inches
print(p)
dev.off()
image_number <<- image_number + 1
## Vulcano Plots ####
### DEP package #### 
# Plot a volcano plot for the contrast of each comparison

for(i in 1:length(comparisons)){
  vulcano<- plot_volcano(dep_analysis, contrast = comparisons[i], label_size = 3, add_names = TRUE, adjusted = FALSE)
  plot(vulcano)
  
  filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_vulcano_DEP_",comparisons[i])  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
  ggsave(paste0(filename,pdf_extension), vulcano, width = 8, height = 6, units = "in")  # Ajusta el tamaño y las unidades según tus necesidades
  ggsave(paste0(filename,tiff_extension), vulcano, width = 8, height = 6, units = "in")  # Ajusta el tamaño y las unidades según tus necesidades
  image_number <<- image_number + 1
}


### ggplot p-value ####
# Define the specific genes to always label

mycolors <- c("blue", "red", "grey")
names(mycolors) <- c("DOWN", "UP", "NO")
for (i in 1:length(comparisons)){
  
  diff_ <- paste(comparisons[i],"diffexpressed",  sep="_")
  ratio <-paste(comparisons[i],"ratio", sep="_")
  pval <- paste(comparisons[i],"p.val", sep="_")
  label <- paste(comparisons[i],"dif_label",sep="_") 
  

  
  vulcanogg <-suppressWarnings(ggplot2::ggplot(data=data_results, 
                                               aes(x=.data[[ratio]],
                                                   y=-log10(.data[[pval]]),
                                                   col=.data[[diff_]])) + 
                                 geom_point() + 
                                 theme_DEP1() + 
                                 ggtitle(paste("Vulcano ", comparisons[i])) +
                                 labs(x=(paste("Fold-Change ",comparisons[i])),
                                      y=expression(paste("-log"[10],"(p-value)"))) +
                                 geom_vline(xintercept= 0) + 
                                 #geom_vline(xintercept= c(-FC, FC), col="red",linetype="dashed") + 
                                 #geom_hline(yintercept= -log10(p_val), col="red",linetype="dashed") + 
                                 scale_colour_manual(values = mycolors) +
                                 geom_text_repel(
                                   aes(label = .data[[label]]),
                                   max.overlaps = 15  # Limit the number of overlaps for additional labels
                                 ) +
                                 guides(colour = FALSE)  # Remove the legend for color
  ) 
  
  print(vulcanogg)
  
  filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_vulcano_ggplot_",comparisons[i])  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
  ggsave(paste0(filename,pdf_extension), vulcanogg, width = 8, height = 6, units = "in")  # Vectorial format
  ggsave(paste0(filename,tiff_extension), vulcanogg, width = 8, height = 6, units = "in")  # Tiff format
  image_number <<- image_number + 1
}

## Barplots ####

  if(print_barplots){

### Adjusted p_value Significant ####


number <- 16# number of graphs per plot.
expressed <- sig_adjusted_data$name


for (j in 1:(length(expressed)/ number)){
  start <- (j-1) * number+1
  end <- j * number
  
  group <- expressed[start:end]
  singleplot <-plot_single(dep_analysis, proteins = group, type="centered") #+ ggtitle(comparisons[i])
  p<-print(singleplot)
  print(p)
  
  filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_enriched_protein")  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
  ggsave(paste0(filename,pdf_extension), p, width = 8, height = 6, units = "in")  # Vectorial format
  ggsave(paste0(filename,tiff_extension), p, width = 8, height = 6, units = "in")  # Tiff format
  image_number <<- image_number + 1
}
### p_value Significant ####
#image_number <- 23
image_number <- image_number+j

number <- 16# number of graphs per plot.
expressed <- significative_data$name

for (j in 1:(length(expressed)/ number)){
  start <- (j-1) * number+1
  end <- j * number
  
  group <- expressed[start:end]
  singleplot <-plot_single(dep_analysis, proteins = group, type="centered") #+ ggtitle(comparisons[i])
  p<-print(singleplot)
  print(p)
  
  filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_enriched_protein_pvalue")  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
  ggsave(paste0(filename,pdf_extension), p, width = 8, height = 8, units = "in")  # Vectorial format
  ggsave(paste0(filename,tiff_extension), p, width = 8, height = 8, units = "in")  # Tiff format
  image_number <<- image_number + 1
}
  }


## Fold-Change correlations ####
ratio_columns <- data_results %>% 
  dplyr::select(ends_with("_ratio")) %>% 
  names()

if(length(comparisons)==1){
  print("No combinations possible. Only one comparison")
}else{
combinations <- combn(ratio_columns, 2)
# Generar gráfico por cada combinación
for (i in 1:ncol(combinations)) {
  x <- combinations[1, i]
  y <- combinations[2, i]
  
  # Crear el gráfico
  plot_title <- paste(ratio_columns[which(ratio_columns == x)], "vs", ratio_columns[which(ratio_columns == y)])
  p <- ggplot(data_results, aes_string(x = x, y = y)) +
    geom_point(shape = 21, fill = "white", size = 2) +
    geom_text_repel(max.overlaps = 15, aes(label = name)) + 
    #geom_mark_hull(aes(fill=affinity),concavity = 2.8)+ #Refine it. Draw only the groups that are being drawn
    labs(title = plot_title, x = ratio_columns[which(ratio_columns == x)], y = ratio_columns[which(ratio_columns == y)])  # Convertir los nombres en símbolos
  
  # Mostrar el gráfico
  print(p)
  filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_FC_correlation",x,"_",y)  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
  ggsave(paste0(filename,pdf_extension), p, width = 8, height = 8, units = "in")  # Vectorial format
  ggsave(paste0(filename,tiff_extension), p, width = 8, height = 8, units = "in")  # Tiff format
  image_number <<- image_number + 1
}
}

#image_number <- 35

# HeatMaps ####

## By comparison ####
### HeatMap padj ####
selected_columns <- paste(comparisons, "_ratio", sep = "")

specific_order <- c("name", selected_columns)

heatmap_df_padj <- data_results %>% 
  dplyr::filter(significant == TRUE) %>% 
  dplyr::select(name, ends_with("_ratio")) %>% 
  #dplyr::arrange(desc(!!!syms(selected_columns))) %>%
  dplyr::select(all_of(specific_order))


# Convert the data to a matrix for plotting

heatmap_matrix_padj <- as.matrix(heatmap_df_padj[,-1])
rownames(heatmap_matrix_padj) <- heatmap_df_padj$name

# Define custom color palette
my_color_palette <- colorRampPalette(viridis(12))(100)#define your color scale

# Plot heatmap with specific order of cluster names on x-axis
heatmap_adj <- pheatmap::pheatmap(heatmap_matrix_padj, cutree_rows = 4, 
                                  cluster_cols = FALSE, 
                                  scale = "none", 
                                  annotation_legend = TRUE,
                                  labels_col = comparisons,  # Add cluster names to x-axis
                                  labels_row = heatmap_df_padj$name,  # Add gene names to y-axis
                                  #annotation_names_col = FALSE,
                                  cluster_rows = TRUE,
                                  main = "Fold-change ratio adjusted-pvalue",
                                  display_numbers = TRUE
)
filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_HeatMap_padj")  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
ggsave(paste0(filename,pdf_extension), heatmap_adj, width = 4, height = 18, units = "in")  # Vectorial format
ggsave(paste0(filename,tiff_extension), heatmap_adj, width = 4, height = 18, units = "in")  # Tiff format

image_number <<- image_number + 1




### HeatMap pval ####

heatmap_df <- data_results %>% 
  dplyr::filter(significance == TRUE) %>% 
  dplyr::select(name,ends_with("_ratio")) %>% 
  # dplyr::arrange(desc(!!!syms(selected_columns))) %>%
  dplyr::select(all_of(specific_order))

# Convert the data to a matrix for plotting

heatmap_matrix <- as.matrix(heatmap_df[,-1])
#rownames(heatmap_matrix) <- heatmap_df$name


# Define custom color palette
my_color_palette <- colorRampPalette(viridis(12))(100)#define your color scale

# Plot heatmap with specific order of cluster names on x-axis
# Plot heatmap with specific order of cluster names on x-axis
heatmap_pval <- pheatmap::pheatmap(heatmap_matrix , 
                                   cluster_rows = TRUE, cutree_rows = 6, 
                                   cluster_cols = FALSE, 
                                   scale = "none", 
                                   labels_col = comparisons,  # Add cluster names to x-axis
                                   annotation_legend = TRUE,
                                   main = "Fold-change ratio",
                                   legend = TRUE,
                                   labels_row = heatmap_df$name,
                                   display_numbers = TRUE
)

filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_HeatMap_ratio")  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
ggsave(paste0(filename,pdf_extension), heatmap_pval, width = 4, height = 18, units = "in")  # Vectorial format
ggsave(paste0(filename,tiff_extension), heatmap_pval, width = 4, height = 18, units = "in")  # Tiff format
image_number <<- image_number + 1


### HeatMap pval excluded####
heatmap_df <- data_results %>% 
  dplyr::filter(significance == TRUE) %>% 
  dplyr::select(name,ends_with("_ratio")) %>% 
  #dplyr::arrange(desc(!!!syms(selected_columns))) %>%
  dplyr::select(all_of(specific_order))

# Convert the data to a matrix for plotting

heatmap_matrix <- as.matrix(heatmap_df[,-1])
rownames(heatmap_matrix) <- heatmap_df$name
# Perform hierarchical clustering
hc <- hclust(dist(heatmap_matrix))  # You may need to adjust the distance metric

# Get cluster assignments
clusters <- cutree(hc, k = 6)  # Assuming you want 6 clusters, adjust as needed

# Add cluster assignments as a new column to heatmap_matrix
heatmap_matrix_with_cluster <- cbind(heatmap_matrix, Cluster = clusters)

#Exclude cluster 6, which it has many non informative proteins. Adjust as required necessities

heatmap_matrix_exclude <- heatmap_matrix_with_cluster[heatmap_matrix_with_cluster[, "Cluster"] != 1, ]
heatmap_matrix_exclude <- heatmap_matrix_exclude[heatmap_matrix_exclude[, "Cluster"] != 5, ]
row_names_ratio <- rownames(heatmap_matrix_exclude)

selected_ratio <- data_results %>% 
  dplyr::filter(name %in% row_names_ratio)

selected_ratio %>% 
  write_xlsx(paste0(output_path,"tables/selected_ratio.xlsx"))


#Remove Cluster Column
cluster_col_index <- which(colnames(heatmap_matrix_exclude) == "Cluster")
heatmap_matrix_exclude <- heatmap_matrix_exclude[, -cluster_col_index]



# Define custom color palette
my_color_palette <- colorRampPalette(viridis(12))(100)#define your color scale
# Plot heatmap with specific order of cluster names on x-axis
heatmap_pval <- pheatmap::pheatmap(heatmap_matrix_exclude, 
                                   cluster_rows = TRUE, cutree_rows = 6,
                                   cluster_cols = FALSE,
                                   scale = "none",
                                   labels_col = comparisons,  # Add cluster names to x-axis
                                   annotation_legend = TRUE,
                                   main = "Fold-change ratio",
                                   display_numbers = TRUE,
                                   labels_row = heatmap_df$name  # Add gene names to y-axis
)

filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_HeatMap_ratio_excluded")  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
ggsave(paste0(filename,pdf_extension), heatmap_pval, width = 4, height = 18, units = "in")  # Vectorial format
ggsave(paste0(filename,tiff_extension), heatmap_pval, width = 4, height = 18, units = "in")  # Tiff format
image_number <<- image_number + 1
#image_number <- 10
## By Log2 Intensity ####

### Log2 Intensity significant####
heatmap_intensity <- data_results %>% 
  dplyr::filter(significant == TRUE) %>% 
  dplyr::select(name,starts_with("lfq_")) 


# Convert the data to a matrix for plotting

heatmap_matrix_intensity <- as.matrix(heatmap_intensity[,-1])
rownames(heatmap_matrix_intensity) <- heatmap_intensity$name

# Define custom color palette
my_color_palette <- colorRampPalette(viridis(12))(100)#define your color scale


# Plot heatmap with specific order of cluster names on x-axis
heatmap_int <- pheatmap::pheatmap(heatmap_matrix_intensity, 
                                  cutree_rows = 4, 
                                  row_names = FALSE, 
                                  scale = "none", 
                                  cluster_cols = FALSE,
                                  #annotation_col = label,
                                  #annotation_colors = list(Gene_Present = c("No" = "white", "Yes" = "black")),
                                  annotation_legend = TRUE,
                                  labels_col = label_columns,  # Add cluster names to x-axis
                                  #labels_row = heatmap_intensity$name,  # Add gene names to y-axis
                                  annotation_names_col = FALSE,
                                  cluster_rows = TRUE,
                                  main = "Log2 Intensity",
                                  legend_title = "Log2 Intensity",  # Add legend title
                                  #color = my_color_palette
)


filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_HeatMap_intensity_padj")  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
ggsave(paste0(filename,pdf_extension), heatmap_int, width = 6, height = 18, units = "in")  # Vectorial format
ggsave(paste0(filename,tiff_extension), heatmap_int, width = 6, height = 18, units = "in")  # Tiff format
image_number <<- image_number + 1




### Log2 Intensity ####
heatmap_intensity <- data_results %>% 
  dplyr::filter(significance == TRUE) %>% 
  dplyr::select(name,starts_with("lfq_")) 


# Convert the data to a matrix for plotting

heatmap_matrix_intensity <- as.matrix(heatmap_intensity[,-1])
#rownames(heatmap_matrix_intensity) <- heatmap_intensity$name

# Define custom color palette
my_color_palette <- colorRampPalette(viridis(12))(100)#define your color scale



# Plot heatmap with specific order of cluster names on x-axis
heatmap_int <- pheatmap::pheatmap(heatmap_matrix_intensity, cutree_rows = 6, row_names = TRUE, 
                                  cluster_cols = FALSE, scale = "none", 
                                  #annotation_col = label,
                                  #annotation_colors = list(Gene_Present = c("No" = "white", "Yes" = "black")),
                                  annotation_legend = TRUE,
                                  labels_col = label_columns,  # Add cluster names to x-axis
                                  #labels_row = heatmap_intensity$name,  # Add gene names to y-axis
                                  annotation_names_col = FALSE,
                                  cluster_rows = TRUE,
                                  main = "Log2 Intensity",
                                  legend_title = "Log2 Intensity",  # Add legend title
                                  #color = my_color_palette,
                                  #labels_row = heatmap_intensity$name,
)


filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_HeatMap_intensity")  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
ggsave(paste0(filename,pdf_extension), heatmap_int, width = 6, height = 18, units = "in")  # Vectorial format
ggsave(paste0(filename,tiff_extension), heatmap_int, width = 6, height = 18, units = "in")  # Tiff format
image_number <<- image_number + 1




### Log2 Intensity Excluded ####
heatmap_intensity <- data_results %>% 
  dplyr::filter(significance == TRUE) %>% 
  dplyr::select(name,starts_with("lfq_")) 



# Convert the data to a matrix for plotting

heatmap_matrix_int_excl <- as.matrix(heatmap_intensity[,-1])
rownames(heatmap_matrix_int_excl) <- heatmap_intensity$name
# Perform hierarchical clustering
hc <- hclust(dist(heatmap_matrix_int_excl))  # You may need to adjust the distance metric

# Get cluster assignments
clusters <- cutree(hc, k = 6)  # Assuming you want 6 clusters, adjust as needed

# Add cluster assignments as a new column to heatmap_matrix
heatmap_matrix_with_cluster_int <- cbind(heatmap_matrix_int_excl, Cluster = clusters)

#Exclude clusters 1&2, which it has many non informative proteins. Adjust as required necessities

heatmap_matrix_exclude_int_1 <- heatmap_matrix_with_cluster_int[heatmap_matrix_with_cluster_int[, "Cluster"] != 1, ]
heatmap_matrix_exclude_int_2 <- heatmap_matrix_exclude_int_1[heatmap_matrix_exclude_int_1[, "Cluster"] != 2, ]

row_names <- rownames(heatmap_matrix_exclude_int_2)

selected_intensity <- data_results %>% 
  dplyr::filter(name %in% row_names)

selected_intensity %>% 
  write_xlsx(paste0(output_path,"tables/selected_intensity.xlsx"))

#Remove Cluster Column
cluster_col_index <- which(colnames(heatmap_matrix_exclude_int_2) == "Cluster")
heatmap_matrix_exclude_int_2 <- heatmap_matrix_exclude_int_2[, -cluster_col_index]



# Define custom color palette
my_color_palette <- colorRampPalette(viridis(12))(100)#define your color scale



# Plot heatmap with specific order of cluster names on x-axis
heatmap_int <- pheatmap::pheatmap(heatmap_matrix_exclude_int_2, cutree_rows = 4, row_names = TRUE, 
                                  cluster_cols = FALSE, scale = "none", 
                                  #annotation_col = label,
                                  #annotation_colors = list(Gene_Present = c("No" = "white", "Yes" = "black")),
                                  annotation_legend = TRUE,
                                  labels_col = label_columns,  # Add cluster names to x-axis
                                  #labels_row = heatmap_matrix_exclude_int,  # Add gene names to y-axis
                                  annotation_names_col = FALSE,
                                  cluster_rows = TRUE,
                                  main = "Log2 Intensity",
                                  legend_title = "Log2 Intensity",  # Add legend title
                                  #color = my_color_palette
)


filename <- paste0(output_path,"figures/", sprintf("%02d", image_number), "_HeatMap_intensity_excluded")  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
ggsave(paste0(filename,pdf_extension), heatmap_int, width = 6, height = 10, units = "in")  # Vectorial format
ggsave(paste0(filename,tiff_extension), heatmap_int, width = 6, height = 10, units = "in")  # Tiff format
image_number <<- image_number + 1



## 





