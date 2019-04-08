# det_opt_cluster_num.R plots the cost function of the k-means algorithm against 
#   the number of clusters. The aim of this procedure it to identify an optimal 
#   number of clusters for k-means clustering. A natural cluster number can be 
#   set to the "elbow" of the plot produced here.

# Author: Ashley Acevedo

library(cluster)
library(factoextra)
source("../cluster_function.R")

# Read table of sequence windows and features:
#   [1] Protein name
#   [2] Start position of sequence window
#   [3] Sequence window 
#   [4 - 23] Fraction of each amino acid (in alphabetical order)
#   [5] Shannon entropy
#   [5:] repeat for exchange groups and higher order ngrams
data = read.table("../processed_sequences/FUS_family_window_50_step_10_ngram_3.txt", 
                  header = TRUE, 
                  stringsAsFactors = FALSE)

ngram.length = 3

# Feature scaling and mean normalization of Entropy
for (n in 1:ngram.length){
  # for amino acid ngrams
  entropy.col = data[[paste("Entropy_", n, sep = "")]]
  data[[paste("Entropy_", n, sep = "")]] = (entropy.col - mean(entropy.col))/
    (max(entropy.col) - min(entropy.col))
}

# Remove columns with no representatives
data.temp = data[ , 4:ncol(data)]
data.temp = data.temp[ , -(which(colSums(data.temp) == 0))]
data = cbind(data[ , 1:3], data.temp)

# For reproducibility
set.seed(1)

# Compute cost and silhouette for a range of cluster numbers
cost = c()
mean.sil = c()
for (i in 2:20){
  kmeans.out = kmeans(data[ , 4:24], centers = i, iter.max = 100, nstart = 100)
  cost = c(cost, kmeans.out$tot.withinss)
  mean.sil = c(mean.sil, mean(silhouette(kmeans.out$cluster, dist(data[ , 4:24]))[ , 3]))
}

# Plot cost versus then number of clusters ("elbow" method)
pdf(file = "Elbow_method.pdf", height = 5, width = 5)
plot(2:20, cost, type = "b", 
     pch  = 16, 
     ylim = c(10, 50),
     xlim = c(2, 20),
     xlab = "K (number of clusters)",
     ylab = "Cost function (sum of square errors)",
     axes = FALSE)
axis(1, at = seq(2, 20, 2), labels = seq(2, 20, 2))
axis(2, at = seq(10, 50, 10), labels = seq(10, 50, 10), las = 1)
dev.off()

# Plot average silhouette versus then number of clusters ("silhouette" method)
pdf(file = "Silhouette_method.pdf", height = 5, width = 5)
plot(2:20, mean.sil, type = "b", 
     pch  = 16, 
     ylim = c(0.1, 0.45),
     xlim = c(2, 20),
     xlab = "K (number of clusters)",
     ylab = "Average silhouette",
     axes = FALSE)
axis(1, at = seq(2, 20, 2), labels = seq(2, 20, 2))
axis(2, at = seq(0.1, 0.45, 0.05), labels = seq(0.1, 0.45, 0.05), las = 1)
dev.off()

# Compute gap statistic
gap.stat = clusGap(data[ , 4:24], FUN = kmeans, nstart = 25, K.max = 20, B = 50)

# Plot gap statistic
pdf(file = "Gap_statistic_method.pdf", height = 5, width = 5)
fviz_gap_stat(gap.stat)
dev.off()



