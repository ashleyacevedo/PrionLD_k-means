# map_clusters_to_gene.R visualizes the mappings of clusters to genes.

# Author: Ashley Acevedo

# Number of clusters
K = 6

# Read table of sequence windows and features:
#   [1] Protein name
#   [2] Start position of sequence window
#   [3] Sequence window 
#   [4 - 23] Fraction of each amino acid (in alphabetical order)
#   [5] Shannon entropy
data = read.table("../processed_sequences/FUS_family_window_50_step_10.txt", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)

# Cluster input data using k-means algorithm
cluster.out = clusterData(data, K, 3, FALSE)
kmeans.out = cluster.out[[1]]
data = cluster.out[[2]]

# Set color scheme for clusters
if (K == 6){
  colors = c("cyan","darkred","orange", "yellow","red", "grey")
} else if (K == 7){
  colors = c("cyan","red","orange", "grey", "darkorange", "yellow","darkred")
} else if (K == 8){
  colors = c("cyan","red","gold","grey","orange", "yellow","darkorange","darkred")
} else { colors = rainbow(K) }

# Map clusters to genes
pdf(file = paste("FUS_family_gene_segment_clusters_K", K, ".pdf", sep = ""), height = 10, width = 7)
par(mar = c(5, 6, 2, 2) )
plot(0, 0, type = "n",xlim = c(0,650), ylim = c(0,22), axes = FALSE, xlab = "", ylab = "")
gene.names = unique(data$ProteinName)
for (n in gene.names){
  c = kmeans.out$cluster[data$ProteinName == n]
  x = data$StartPosition[data$ProteinName == n]
  for (i in 1:length(x)){
    rect(x[i], 21 - which(gene.names == n) - 0.25, x[i] + nchar(data$Sequence[1]), 21 - which(gene.names == n) + 0.25,
         col = colors[c[i]], border = NA)
  }
}

gene.names = strsplit(gene.names, split = "_")
gene.names = unlist(gene.names)
gene.names = gene.names[seq(1, length(gene.names) - 1, 2)]

axis(1, at = seq(0, 600, 100), labels = seq(0, 600, 100))
axis(2, at = 0:20, labels = rev(gene.names), tick = FALSE, las = 1)
legend(500, 5, legend = as.character(1:K), col = colors, bty = "n", title = "cluster", pch = 15)
dev.off()