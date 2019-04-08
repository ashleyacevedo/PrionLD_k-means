# pca_cluster_viz.R reduces dimensionality of sequence features so that
#   clusters can be visualized in 2D.

# Author: Ashley Acevedo

# Load required packages
library(ggplot2)
source("../cluster_function.R")

# Number of clusters
K = 6

# Read table of sequence windows and features:
#   [1] Protein name
#   [2] Start position of sequence window
#   [3] Sequence window 
#   [4:23] Fraction of each amino acid (in alphabetical order)
#   [5] Shannon entropy
#   [5:] repeat for exchange groups and higher order ngrams
data = read.table("../processed_sequences/FUS_family_window_50_step_10_ngram_3.txt", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)

# Cluster input data using k-means algorithm
cluster.out = clusterData(data, K, 3, FALSE)
kmeans.out = cluster.out[[1]]
data = cluster.out[[2]]

# Reduce dimensions to support 2D visualization
pca.out = prcomp(data[ , 4:ncol(data)], center = TRUE, scale = TRUE)
percent.exp.var <- round(pca.out$sdev^2/sum(pca.out$sdev^2)*100)
pca.loadings = data.frame(variables = rownames(pca.out$rotation), pca.out$rotation)

# Visualize clusters 
pca.df = data.frame("PC1" = pca.out$x[ , 1],
                    "PC2" = pca.out$x[ , 2],
                    "PC3" = pca.out$x[ , 3],
                    "cluster" = as.factor(kmeans.out$cluster))

# Set color scheme for clusters
if (K == 6){
  colors = c("cyan","darkred","orange", "yellow","red", "grey")
} else if (K == 7){
  colors = c("cyan","red","orange", "grey", "darkorange", "yellow","darkred")
} else if (K == 8){
  colors = c("cyan","red","gold","grey","orange", "yellow","darkorange","darkred")
} else { colors = rainbow(K) }


pdf(file = paste("PCA_cluster_visualization_PC1_vs_PC2_K", K, ".pdf", sep = ""), height = 5, width = 5.5)
ggplot(pca.df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point() + 
  xlab(paste("PC1 (", percent.exp.var[1], "%)")) +
  ylab(paste("PC2 (", percent.exp.var[2], "%)")) +
  scale_color_manual(values = colors)
#  geom_segment(data = pca.loadings, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), 
#               arrow = arrow(length = unit(1/2, "picas")),
#               color = "darkgrey", size = 0.75) +
#  annotate("text", x = (pca.loadings$PC1*8.5), y = (pca.loadings$PC2*8.5),
#           label = pca.loadings$variables)
dev.off()

pdf(file = paste("PCA_cluster_visualization_PC2_vs_PC3_K", K, ".pdf", sep = ""), height = 5, width = 5.5)
ggplot(pca.df, aes(x = PC3, y = PC2, color = cluster)) +
  geom_point() + 
  xlab(paste("PC3 (", percent.exp.var[3], "%)")) +
  ylab(paste("PC2 (", percent.exp.var[2], "%)")) +
  scale_color_manual(values = colors)
#  geom_segment(data = pca.loadings, aes(x = 0, y = 0, xend = (PC3*8), yend = (PC2*8)), 
#               arrow = arrow(length = unit(1/2, "picas")),
#               color = "darkgrey", size = 0.75) +
#  annotate("text", x = (pca.loadings$PC3*8.5), y = (pca.loadings$PC2*8.5),
#           label = pca.loadings$variables)
dev.off()

pdf(file = paste("PCA_cluster_visualization_PC1_vs_PC3_K", K, ".pdf", sep = ""), height = 5, width = 5.5)
ggplot(pca.df, aes(x = PC3, y = PC1, color = cluster)) +
  geom_point() + 
  xlab(paste("PC3 (", percent.exp.var[3], "%)")) +
  ylab(paste("PC1 (", percent.exp.var[1], "%)")) +
  scale_color_manual(values = colors)
#  geom_segment(data = pca.loadings, aes(x = 0, y = 0, xend = (PC3*8), yend = (PC1*8)), 
#               arrow = arrow(length = unit(1/2, "picas")),
#               color = "darkgrey", size = 0.75) +
#  annotate("text", x = (pca.loadings$PC3*8.5), y = (pca.loadings$PC1*8.5),
#           label = pca.loadings$variables)
dev.off()

