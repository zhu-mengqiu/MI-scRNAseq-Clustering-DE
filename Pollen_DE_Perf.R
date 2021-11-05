create_de_plot <- function(ds_list, name, pow, nimpute){
  p.im = ds_list[[name]]
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

path.pol <- '~/SingleCell/Pollen/Pollen_'
path.kol <- '~/SingleCell/Kolo/Kolo_'
path.dar <- '~/SingleCell/Darmanis/Darmanis_MI/Darmanis_'

nimpute <- 100
p.im.pol <- readRDS(paste0(path.pol, as.character(nimpute), '_kw_im_p.RDS'))
p.im.kol <- readRDS(paste0(path.kol, as.character(nimpute), '_kw_im_p.RDS'))
p.im.dar <- readRDS(paste0(path.dar, as.character(nimpute), '_kw_im_p.RDS'))

ds_list <- list('Pollen' = p.im.pol, 
                'Kolodziejczyk' = p.im.kol, 
                'Darmanis' = p.im.dar)

create_de_plot(ds_list, 'Pollen', 1/6, nimpute)
create_de_plot(ds_list, 'Kolodziejczyk', 1/4, nimpute)
create_de_plot(ds_list, 'Darmanis', 1/3, nimpute)
