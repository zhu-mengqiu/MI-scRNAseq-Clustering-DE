# set working directory to the folder where output data from Pollen_Clust.R 
# is stored

# name of output data is Pollen_n_KW_P.RDS
# n = number of imputations

#----------FDR Plot----------
# define a function to create FDR plot
create_de_plot <- function(ds, name, pow, nimpute){
  p.im = ds
  p.mi <- apply(p.im, 1, median)
  
  p.im.adj <- apply(p.im, 2, function(c) p.adjust(c, method = 'BH'))
  p.mi.adj <- p.adjust(p.mi, method = 'BH')
  
  p.im.adj.ord <- apply(p.im.adj, 2, function(c) c[order(c)])
  p.mi.adj.ord <- p.mi.adj[order(p.mi.adj)]
  
  len <- length(p.mi.adj.ord)
  seq.x <- seq(1, len %/% 5000) * 5000
  x <- seq(1,len,1)
  x.1 <- x^6
  
  seq.y <- seq(0, 1, 0.2)
  seq.y.1 <- seq.y^pow
  
  setEPS()
  postscript(paste0('DE_plot_',name,'_', as.character(nimpute),'.eps'), 
             width=12, height=12)
  
  par(oma=c(1,1,1,1))
  par(mar=c(5,5,4,2))
  
  y <- p.im.adj.ord[,1]
  y.1 <- y^pow
  
  plot(x.1, 
       y.1, 
       type = 'l', 
       col = 'grey',
       ylim = c(0,1),
       ylab = 'FDR',
       xlab = 'Number of Differentially Expressed Genes',
       xaxt = 'n',
       yaxt = 'n',
       lty = 1,
       lwd = 1,
       main = name,
       cex.main = 2.5,
       cex.lab = 2)
  
  axis(1, at=x.1[seq.x], labels=x[seq.x], cex.axis=1.5)
  axis(2, at=c(seq.y.1, 0.05^pow, 0.01^pow, 0.1^pow), 
       labels=c(seq.y, 0.05, 0.01, 0.1), cex.axis=1.5)
  
  for (i in 2:nimpute){
    y <- p.im.adj.ord[,i]
    y.1 <- y^pow
    
    lines(x.1,
          y.1, 
          type = 'l', 
          col = 'grey',
          lty = 1, 
          lwd = 1)
  }
  
  y <- p.mi.adj.ord
  y.1 <- y^pow
  
  lines(x.1,
        y.1,
        type = 'l',
        col = 'black',
        lty = 4,
        lwd = 5)
  
  abline(h=0.05^pow, col='black', lty=2, lwd=3)
  abline(h=0.01^pow, col='black', lty=2, lwd=3)
  abline(h=0.1^pow, col='black', lty=2, lwd=3)
  
  legend('topleft',
         legend = c('Multiple Imputation', 'One-time Imputation'),
         bty = 'n',
         col = c('black', 'grey'),
         cex = 2,
         lty = c(4,1),
         lwd = c(5,5))
  
  dev.off()
}

nimpute <- 100
p.im.pol <- readRDS(paste0('Pollen_',as.character(nimpute), '_KW_P.RDS'))

create_de_plot(p.im.pol, 'Pollen', 1/6, nimpute)

#----------MI Rank----------
# define a function to compute rank of MI results
compute_mi_rank <- function(p.im, p.mi, alpha){
  p.im.adj <- apply(p.im, 2, function(c) p.adjust(c, method = 'BH'))
  p.mi.adj <- p.adjust(p.mi, method = 'BH')
  rej.im <- apply(p.im.adj, 2, function(c) sum(c <= alpha))
  rej.mi <- sum(p.mi.adj <= alpha)
  rank <- sum(rej.mi < rej.im) + 1
  return (rank)
}

p.mi.pol <- apply(p.im.pol, 1, median)
print(compute_mi_rank(p.im.pol, p.mi.pol, 0.01))
print(compute_mi_rank(p.im.pol, p.mi.pol, 0.05))
print(compute_mi_rank(p.im.pol, p.mi.pol, 0.1))
