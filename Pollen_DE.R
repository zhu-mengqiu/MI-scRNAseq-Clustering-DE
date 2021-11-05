# set working directory to the folder where imputed datasets are stored

# names of imputed datasets set as Pollen_scIGANs_impute_nolabel_n.txt
# n = 1,2,3....

# set seed for reproducible results
set.seed(200)

# set number of imputations for multiple imputation-based clustering
nimpute <- 100

# cell type labels assigned by authors
true.label <- read.table('Pollen_label.txt')$V1
ncluster <- length(unique(true.label))

# define a function to preprocess data
# normalization: divide each cell by sum of counts and multiply by one million
# log transformation: take log of (count + 1)
preprocess <- function(m){
  m <- sweep(m, 2, colSums(m), '/') * 1e6
  m <- log10(m+1)
  return(m)
}

# define a function to return p-value of Kruskal-Wallis test
kwp <- function(count, label){
  test <- kruskal.test(count ~ as.factor(label))
  p <- unname(test$p.value)
  return (p)
}

#----------Imputed Data----------
data.im <- vector(mode = 'list', length = nimpute)
for (i in 1:nimpute){
  data.file <- paste0('Pollen_scIGANs_impute_nolabel_', as.character(i), 
                      '.txt')
  data.tmp <- read.table(data.file, header = T, sep = '\t')
  data.im[[i]] <- preprocess(as.matrix(data.tmp[,-1]))
}

#----------Individual DE Kruskal-Wallis----------
ngene <- dim(data.im[[1]])[1]
kw.im.p <- matrix(data = NA, nrow = ngene, ncol = nimpute)
for (i in 1:nimpute){
  kw.im.p[,i] <- apply(data.im[[i]], 1, function(r) kwp(r, true.label))
}

#----------Output----------
saveRDS(object = kw.im.p, file = paste0('Pollen_', as.character(nimpute), 
                                        '_KW_P.RDS'))