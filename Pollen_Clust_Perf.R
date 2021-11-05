# set working directory to the folder where output data from Pollen_Clust.R 
# is stored

# name of output data is Pollen_n_KMeans_Louvain.RDS
# n = number of imputations

# load packages
library(xtable)

# number of imputations
nimpute <- 100

# performance metrics
nmetric <- 7
metrics <- c('NMI', 'ARI', 'NVI', 'NID', 'ACC', 'AUC', 'F.score')

# load performance metric values
tab.pol <- readRDS('Pollen_100_KMeans_Louvain.RDS')

#----------Performance Metrics: Table----------
# define a function to create a table of performance metrics
perf.table <- function(tab){
  tab <- tab[,metrics]
  SI <- tab[-(nimpute+1),]
  MI <- tab[(nimpute+1),]
  SI.Quantile <- apply(tab[-(nimpute+1),], 2, quantile)[2:4,]
  MI.Percentile <- apply(tab, 2, 
                         function(c) mean(c[-(nimpute+1)] < c[(nimpute+1)])*100)
  perf <- rbind(SI.Quantile, MI, MI.Percentile)
  rownames(perf) <- c('Q1', 'Median', 'Q3', 'MI', 'Percent')
  perf.tab <- xtable(perf, digits = 3)
  return (perf.tab)
}

# print the table
perf.pol <- perf.table(tab.pol)
print(perf.pol)

#----------Performance Metrics: Histograms----------
setEPS()
postscript(paste0('Clust_Perf_Hist.eps'), width=5, height=21)
par(mfcol=c(7,1))
par(oma=c(1,1,1,1))
par(mar=c(5,5,4,2))
tab.pol <- tab.pol[,metrics]
for (metric.name in metrics){
  metric <- tab.pol[,metric.name]
  metric.si <- metric[-(nimpute+1)]
  metric.mi <- metric[nimpute+1]
  hist(metric.si,
       main=paste0('Pollen_', metric.name),
       breaks=10,
       xlab=metric.name,
       cex.main=2,
       cex.lab=2,
       font.main=2,
       font.lab=2)
  abline(v=metric.mi, lty=2, lwd=4)
}
dev.off()

#----------Performance Metrics: ECDF Plots----------
setEPS()
postscript(paste0('Clust_Perf_ECDF.eps'), width=5, height=21)
par(mfcol=c(7,1))
par(oma=c(1,1,1,1))
par(mar=c(5,5,4,2))
tab.pol <- tab.pol[,metrics]
for (metric.name in metrics){
  metric <- tab.pol[,metric.name]
  metric.si <- metric[-(nimpute+1)]
  metric.mi <- metric[nimpute+1]
  ecdf.si <- ecdf(metric.si)
  plot(ecdf.si,
       main=paste0('Pollen_', metric.name),
       xlab=metric.name,
       cex.main=2,
       cex.lab=2,
       font.main=2,
       font.lab=2)
  abline(v=metric.mi, lty=2, lwd=4)
}
dev.off()