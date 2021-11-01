setwd('/CCAS/home/mzhu32/Pollen')
# setwd('~/SingleCell/Pollen')

set.seed(200)
nimpute <- 20

# true label
true.label <- read.table('Pollen_label.txt')$V1
ncluster <- length(unique(true.label))

# define a function to preprocess data
preprocess <- function(m){
  m <- sweep(m, 2, colSums(m), '/') * 1e6
  m <- log10(m+1)
  return(m)
}

# define a function to return F and p-value of anova test
aovt <- function(count, label){
  test <- anova(aov(count ~ as.factor(label)))
  s <- test$`F value`[1]
  p <- test$`Pr(>F)`[1]
  return (list('s'=s, 'p'=p))
}

# define a function to return F and p-value of kw test
kwt <- function(count, label){
  test <- kruskal.test(count ~ as.factor(label))
  s <- unname(test$statistic)
  p <- unname(test$p.value)
  return (list('s'=s, 'p'=p))
}

#----------raw data----------
data.raw <- read.table('Pollen_raw.txt', header = T, sep = '\t')
genes <- data.raw$GENE_ID
data.raw <- as.matrix(data.raw[,-1])
ngene <- dim(data.raw)[1]
ncell <- dim(data.raw)[2]
data.raw <- preprocess(data.raw)

#----------imputed data----------
data.im <- vector(mode = 'list', length = nimpute)
for (i in 1:nimpute){
  data.file <- paste0('Pollen_scIGANs_impute_nolabel_', as.character(i), 
                      '.txt')
  data.tmp <- read.table(data.file, header = T, sep = '\t')
  data.im[[i]] <- preprocess(as.matrix(data.tmp[,-1]))
}

#----------DE-Raw----------
# aovsp.raw <- matrix(
#   unlist(apply(data.raw, 1, function(r) aovt(r, true.label))),
#   ncol = 2,
#   byrow = T)

# kwsp.raw <- matrix(
#   unlist(apply(data.raw, 1, function(r) kwt(r, true.label))),
#   ncol = 2,
#   byrow = T)
  
#----------DE-Individual----------
# aovsp.im.s <- matrix(data = NA, nrow = ngene, ncol = nimpute)
# aovsp.im.p <- matrix(data = NA, nrow = ngene, ncol = nimpute)
# for (i in 1:nimpute){
#   aovsp.im <- matrix(
#     unlist(apply(data.im[[i]], 1, function(r) aovt(r, true.label))),
#     ncol = 2,
#     byrow = T)
#   aovsp.im.s[,i] <- aovsp.im[,1]
#   aovsp.im.p[,i] <- aovsp.im[,2]
# }

kwsp.im.s <- matrix(data = NA, nrow = ngene, ncol = nimpute)
kwsp.im.p <- matrix(data = NA, nrow = ngene, ncol = nimpute)
for (i in 1:nimpute){
  kwsp.im <- matrix(
    unlist(apply(data.im[[i]], 1, function(r) kwt(r, true.label))),
    ncol = 2,
    byrow = T)
  kwsp.im.s[,i] <- kwsp.im[,1]
  kwsp.im.p[,i] <- kwsp.im[,2]
}

#----------Output----------
# saveRDS(object = aovsp.raw[,1], file = 'Pollen_aov_raw_s.RDS')
# saveRDS(object = aovsp.raw[,2], file = 'Pollen_aov_raw_p.RDS')
# saveRDS(object = aovsp.im.s, file = 'Pollen_aov_im_s.RDS')
# saveRDS(object = aovsp.im.p, file = 'Pollen_aov_im_p.RDS')

# saveRDS(object = kwsp.raw[,1], file = 'Pollen_kw_raw_s.RDS')
# saveRDS(object = kwsp.raw[,2], file = 'Pollen_kw_raw_p.RDS')
saveRDS(object = kwsp.im.s, file = paste0('Pollen_', as.character(nimpute), 
                                          '_kw_im_s.RDS'))
saveRDS(object = kwsp.im.p, file = paste0('Pollen_', as.character(nimpute), 
                                          '_kw_im_p.RDS'))