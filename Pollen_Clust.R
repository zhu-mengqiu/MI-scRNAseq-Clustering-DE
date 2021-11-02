setwd('/CCAS/home/mzhu32/Pollen')

library(ROCR)
library(aricode)
library(scran)
library(igraph)

set.seed(200)
nimpute <- 100

# true label
true.label <- read.table('Pollen_label.txt')$V1
ncluster <- length(unique(true.label))

# define a function to preprocess data
# normalization: divide each cell by sum of counts and multiply by 10K
# log transformation: take log of (count + 1)
preprocess <- function(m){
  m <- sweep(m, 2, colSums(m), '/') * 1e6
  m <- log10(m+1)
  return(m)
}

# define a function to compute classification metrics: 
#   accuracy, roc auc, and f score
# positive class: two cells belong to the same cell type
# negative class: two cells belong to different cell types

class.metrics <- function(c1, c2){
  # c1, true labels
  # c2, predicted labels
  ref.label <- as.vector(outer(c1, c1, '=='))
  pred.label <- as.vector(outer(c2, c2, '=='))
  pred <- ROCR::prediction(as.numeric(pred.label), ref.label)
  
  # accuracy
  acc.tmp <- ROCR::performance(pred, "acc")
  acc <- as.numeric(acc.tmp@y.values[[1]][2])
  
  # ROC area under the curve
  auc.tmp <- ROCR::performance(pred, "auc")
  auc <- as.numeric(auc.tmp@y.values)
  
  # F1 score
  f.tmp <- performance(pred, "f")
  f <- as.numeric(f.tmp@y.values[[1]][2])
  
  # output
  return(list(F.score = f, AUC = auc, ACC = acc))
}

# define a function to combine individual clusterings to a consensus matrix
consmx <- function(clusts, method){
  n <- dim(clusts)[2]
  bimx <- vector(mode = 'list', length = n)
  for (i in 1:n){
    label <- clusts[,i]
    bimx[[i]] <- outer(label, label, '==')
  }
  consmx <- apply(simplify2array(bimx), 1:2, method)
  consmx <- apply(consmx, 1:2, as.numeric) 
}

#----------Imputed Data----------
data.im <- vector(mode = 'list', length = nimpute)
for (i in 1:nimpute){
  data.file <- paste0('Pollen_scIGANs_impute_nolabel_', as.character(i), 
                      '.txt')
  data <- read.table(data.file, header = T, sep = '\t')
  data.im[[i]] <- preprocess(as.matrix(data[,-1]))
}

# ----------Individual Clustering K-Means----------
nrun <- 1

clusts <- vector(mode = 'list', length = nimpute)
for (i in 1:nimpute){
  label <- kmeans(t(data.im[[i]]), centers = ncluster, nstart = nrun)$cluster
  clusts[[i]] <- unlist(label)
}
clusts <- do.call(cbind, clusts)

#----------Performance Single Imputation----------
SI <- vector(mode = 'list', length = nimpute)
for (i in 1:nimpute){
  im.label <- clusts[,i]
  perf.class <- class.metrics(true.label, im.label)
  perf.clust <- aricode::clustComp(true.label, im.label)
  SI[[i]] <- unlist(c(perf.clust, perf.class))
}

SI <- do.call(rbind, SI)

#----------Consensus Clustering Louvain----------
method <- 'median'
consensus <- consmx(clusts, method)
consensus.graph <- graph.adjacency(consensus, 
                                   mode = 'undirected', 
                                   weighted = NULL,
                                   diag = F)
mi.label <- igraph::cluster_louvain(consensus.graph)$membership  

#----------Performance Multiple Imputation----------
perf.class <- class.metrics(true.label, mi.label)
perf.clust <- aricode::clustComp(true.label, mi.label)
MI <- unlist(c(perf.clust, perf.class))

#----------Output----------
Perf <- rbind(SI, MI)
saveRDS(Perf, file = paste0('Pollen_', as.character(nimpute),
                            '_noScale_KMeans_', as.character(nrun),
                            '_Louvain_', method,'.RDS'))