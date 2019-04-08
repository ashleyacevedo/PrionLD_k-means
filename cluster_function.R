# cluster_function.R prepares the ngram encoded feature table by scaling and 
#   normalizing the entropy columns and removing columns with no representatives
#   in the table. Once the table is prepared, k-means clustering is performed.
#   The function takes as input, an initial feature table, the number of 
#   desired clusters, the maximum ngram length and a logical for the presence
#   of ngrams based on exchange groups. The function returns the output of k-
#   means clustering and the cleaned feature table.

clusterData = function(data, K, ngram.length, exchange){
  
  # Feature scaling and mean normalization of Entropy
  for (n in 1:ngram.length){
    
    # for amino acid ngrams
    entropy.col = data[[paste("Entropy_", n, sep = "")]]
    data[[paste("Entropy_", n, sep = "")]] = (entropy.col - mean(entropy.col))/
                                             (max(entropy.col) - min(entropy.col))
    
    # for exchange group ngrams
    if (exchange){
    entropy.col = data[[paste("exEntropy_", n, sep = "")]]
    data[[paste("exEntropy_", n, sep = "")]] = (entropy.col - mean(entropy.col))/
                                               (max(entropy.col) - min(entropy.col))
    }
  }
  
  # Remove columns with no representatives
  data.temp = data[ , 4:ncol(data)]
  data.temp = data.temp[ , -(which(colSums(data.temp) == 0))]
  data = cbind(data[ , 1:3], data.temp)
  
  # For reproducibility
  set.seed(1)
  
  # Assign clusters
  kmeans.out = kmeans(data[ , 4:ncol(data)], centers = K, iter.max = 100, nstart = 100) 

  list(kmeans.out, data)
}

