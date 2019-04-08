# label_PLDs.R assigns PLD label to FUS family genes and outputs a new file
#   with sequence features for only those sequences labeled as PLDs. Some genes
#   are ignored because their most likely PLD domain is either much longer than
#   expected or not in the expected location as shown in Wang et al. 2018 Cell.

# Author: Ashley Acevedo

# Load required packages
source("../cluster_function.R")
library(ggplot2)

# Number of clusters
K = 7

# Read table of sequence windows and features:
#   [1] Protein name
#   [2] Start position of sequence window
#   [3] Sequence window 
#   [4:23] Fraction of each amino acid (in alphabetical order)
#   [5] Shannon entropy
#   [5:] repeat for exchange groups and higher order ngrams
data.in = read.table("../processed_sequences/FUS_family_window_50_step_10_ngram_3.txt", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)

# Cluster input data using k-means algorithm
cluster.out = clusterData(data.in, K, 3, FALSE)
kmeans.out = cluster.out[[1]]
data = cluster.out[[2]]

# Genes for which PLDs can be labeled
genes = unique(data$ProteinName)
#   The longest stretch of sequence mapping to suspected pld does not
#   match what is shown in figure 1 of Wang et al. 2018 Cell
genes = genes[-c(which(genes == "hnRNPH3_P31942"), 
                 which(genes == "CSTF2_P33240"),
                 which(genes == "CSTF2T_Q9H0L4"))]

# Find the longest stretch of sequence mapping to a pld cluster for
# each gene. Save the indicies of those sequences.
plds = c()
for (gene in genes){
  inds = which(data$ProteinName == gene)
  clus = kmeans.out$cluster[inds]
  inds = inds[clus != 1 & clus != 4]
  inds.cs = cumsum(c(1, diff(inds) - 1))
  runs = rle(inds.cs)
  inds = inds[which(inds.cs == runs[["values"]][which.max(runs[["lengths"]])])]
  plds = c(plds, inds)
}

label = rep("other", nrow(data))
label[plds] = "pld"

# Reduce dimensions to support 2D visualization
pca.out = prcomp(data[ , 4:ncol(data)], center = TRUE, scale = TRUE)
percent.exp.var <- round(pca.out$sdev^2/sum(pca.out$sdev^2)*100)
pca.loadings = data.frame(variables = rownames(pca.out$rotation), pca.out$rotation)

# Visualize clusters with labels
pca.df = data.frame("PC1" = pca.out$x[ , 1],
                    "PC2" = pca.out$x[ , 2],
                    "PC3" = pca.out$x[ , 3],
                    "label" = as.factor(label))

pdf(file = "PCA_cluster_visualization_PC1_vs_PC2_K7_labeled_seqs - ngram 3.pdf", height = 5, width = 5.5)
ggplot(pca.df, aes(x = PC1, y = PC2, color = label)) +
  geom_point() + 
  xlab(paste("PC1 (", percent.exp.var[1], "%)")) +
  ylab(paste("PC2 (", percent.exp.var[2], "%)")) +
  scale_color_manual(values = c("grey", "dodgerblue"))
dev.off()

# Compute of each sample from its assigned centroid
cent.d = c()
for (i in 1:nrow(data)){
  clus = kmeans.out$cluster[i]
  d = dist(rbind(data[i, 4:ncol(data)], kmeans.out$centers[clus, ]), method = "euclidean")
  cent.d = c(cent.d, d)
}

# Compute the average distance of samples from each centroid for sequences
# assigned as PLDs and sequences not assigned as PLDs. Write to a file.
pld.clus = c(2, 3, 5, 6, 7)
d.other = c()
d.pld = c()
for (clus in pld.clus){
  d.other = c(d.other, round(mean(cent.d[kmeans.out$cluster == clus & label == "other"]), 3))
  d.pld = c(d.pld, round(mean(cent.d[kmeans.out$cluster == clus & label == "pld"]), 3))
}
distance.df = data.frame("cluster" = pld.clus, "distance.non.pld" = d.other, "distance.pld" = d.pld)

# Write table of average distances from each cluster centroid of samples
#   assigned to that centroid
write.table(distance.df,"Distance from centroid - non-PLD v PLD - K7 - ngram 3.txt", 
            sep="\t", row.names=FALSE, quote = FALSE)

# Write new file with features for only the sequences labeled as PLDs
write.table(data.in[plds,], "FUS_family_window_50_step_10_ngram_3_pld_only.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)