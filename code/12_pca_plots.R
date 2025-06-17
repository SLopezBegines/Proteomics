

# libraries ####
source("../code/00_packages.R")
source("../code/global_variables.R")


# Load data ####
#Load dep_analysis object
load(paste0(output_path,"RData/dep_analysis.RData"))
#Load data results tables
data_results <- readxl::read_xlsx(path = paste0(output_path,"tables/data_results.xlsx"))
significative_data <- readxl::read_xlsx(path = paste0(output_path,"tables/significative_data.xlsx"))
sig_adjusted_data <- readxl::read_xlsx(path = paste0(output_path,"tables/sig_adjusted_data.xlsx"))


label <-    Exp_design$columns_to_rename


LFQ_columns <- c(grep(label[1],colnames(data_results))[1]:((grep(label[1],colnames(data_results))[1]+length(label))-1))
# get LFQ column numbers
LFQ_columns <- data_results %>% 
  dplyr::select(starts_with("lfq_")) %>% 
  colnames()


pca_matrix <- data_results[,c("significant", "significance",LFQ_columns)]
pca_matrix <- as.matrix(data_results[,c("significant", "significance",LFQ_columns)])
row.names(pca_matrix) <- data_results$ID

pca_res <- prcomp(pca_matrix, scale. = TRUE)
pca_res$group       <- c(rep("WT", 3), rep("CLN3_Lux1", 3), rep("CLN3_Lux2", 3), rep("CLN12",3))
pca_res$repetition <- rep(c(rep("1", 1), rep("2", 1),rep("3", 1)), 4)
#loadings by adj_p_value
pca <- autoplot(pca_res, data = data_results, colour= 'significant',shape = TRUE)
print(pca)
filename <- paste0(output_path,"figures/",image_number+1,"_PCA_loadings_p_adj")  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
ggsave(paste0(filename,pdf_extension), pca, width = 5, height = 5, units = "in")  # Vectorial format
ggsave(paste0(filename,tiff_extension), pca, width = 5, height = 5, units = "in")  # Tiff format
image_number <- image_number+1
#loadings by p_value
pca <- autoplot(pca_res, data = data_results, colour= 'significance',shape = TRUE)
print(pca)
filename <- paste0(output_path,"figures/",image_number+1,"_PCA_loadings_p_adj")  # Nombre del archivo de salida (puedes cambiar la extensión según el formato deseado)
ggsave(paste0(filename,pdf_extension), pca, width = 5, height = 5, units = "in")  # Vectorial format
ggsave(paste0(filename,tiff_extension), pca, width = 5, height = 5, units = "in")  # Tiff format
image_number <- image_number+1


autoplot(pca_res,data = data_results, colour= 'group',
         loadings = TRUE, loadings.colour = pca_res$group,
         loadings.label = TRUE, loadings.label.size = 3)


pca <- autoplot(kmeans(pca_matrix,4), data = data_results,colour= 'significant',
         loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3)
print(pca)

#########
#https://mda.tools/docs/news.html

pca_matrix <- as.matrix(data_results[,LFQ_columns])
row.names(pca_matrix) <- data_results$ID
m = pca(pca_matrix, 7, scale = TRUE, info = "People PCA model")
m = selectCompNum(m, 5)


par(mfrow = c(1, 2))
mdaplot(m$res$cal$scores, type = "p", show.labels = FALSE, show.lines = c(0, 0))
mdaplot(m$loadings, type = "p", show.labels = TRUE, show.lines = c(0, 0))

#####
#BiocManager::install("edgeR")
library(edgeR)


pca_matrix <- as.matrix(data_results[,LFQ_columns])
row.names(pca_matrix) <- data_results$ID

# Make a DGEList and add metadata
y <- DGEList(counts = pca_matrix)
y$samples$group       <- c(rep("WT", 3), rep("CLN3_Lux1", 3), rep("CLN3_Lux2", 3), rep("CLN12",3))
y$samples$repetition <- rep(c(rep("1", 1), rep("2", 1),rep("3", 1)), 4)

# Normalize and obtain logcounts for QC
y <- calcNormFactors(y)
logCPMs <- cpm(y, log = TRUE)


# Calculate rowwise variance
rv <- apply(logCPMs, 1, var)

# Sort decreasingly and take top 1000
o <- order(rv, decreasing=TRUE)
top1000 <- head(o, 1000)

# From the logCPMs subset for the top-1000
logCPM_top1000 <- logCPMs[top1000,]

# Run PCA
pca <- prcomp(t(logCPM_top1000))

# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, y$samples)

# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100

# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(x=PC1, y=PC2, color=group, shape=repetition)) + 
  geom_point(size=3) +
  xlab(labs[1]) + ylab(labs[2])


# correct for the batch which here is the "kit"
batch <- factor(y$samples$group)

logCPMs_corrected <- limma::removeBatchEffect(logCPMs, batch = batch)

# repeat PCA as before, using the same genes
logCPM_corrected_top1000 <- logCPMs_corrected[top1000,]

# Run PCA
pca <- prcomp(t(logCPM_corrected_top1000))

# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x, y$samples)

# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100

# We focus here on PC1 and PC2
use.pcs <- c(1,2)
labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

ggplot(to_plot, aes(x=PC1, y=PC2, color=group, shape=repetition)) + 
  geom_point(size=3) +
  xlab(labs[1]) + ylab(labs[2])



#####
#Differential expression

# Design accounting for kit (=batch)
design <- model.matrix(~group+repetition, y$samples)

# QLF workflow from edgeR
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# see head(design) => the fourth column is treatment which is what we want to test
fit  <- glmQLFTest(fit, coef = 4)

# get stats as a data.frame
tt <- data.frame(topTags(fit, n=Inf))

# Classify genes into significantly up and down
tt_modified <- tt %>% 
  mutate(status=factor(case_when(logFC>0 & FDR<0.05 ~ "up",
                                 logFC<0 & FDR<0.05 ~ "down",
                                 TRUE ~ "not.signif"),
                       levels=c("up", "not.signif", "down")))

# MA-plot
ggplot(tt_modified, aes(x=logCPM, y=logFC, color=status)) +
  geom_point(size=1) +
  scale_color_manual(values=c("firebrick", "grey", "dodgerblue")) +
  ggtitle("MA plot")

# Volcano (logFC vs -log10(pvalue -- I prefer FDR))
ggplot(tt_modified, aes(x=logFC, y=-log10(FDR), color=status)) +
  geom_point(size=1) +
  scale_color_manual(values=c("firebrick", "grey", "dodgerblue")) +
  ggtitle("Volcano-plot")


###########
pca_matrix <- data_results[,c("significant", "significance",LFQ_columns)]

# Prepare data - Remove non-numeric columns (ID, significant, significance)
data <- pca_matrix[, -(1:3)]

# Load necessary libraries
library(ggplot2)
library(ggrepel)  # For geom_text_repel



# Compute PCA
pca_result <- prcomp(data, scale. = TRUE)

# Create PCA plot for scores
scores_plot <- autoplot(pca_result, data = pca_matrix, colour = 'significant', shape = 'significance', label = TRUE, label.size = 3) +
  theme_minimal() +
  ggtitle("PCA Scores Plot")

# Create PCA plot for loadings
loadings_plot <- autoplot(pca_result, data = pca_matrix, loadings = TRUE, loadings.label = TRUE) +
  theme_minimal() +
  ggtitle("PCA Loadings Plot")

# Display plots
print(scores_plot)
print(loadings_plot)
