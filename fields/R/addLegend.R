addLegend<- function(legendLayout,
                col,
               zlim,
          axis.args = NULL,
        legend.args = NULL,
         legend.cex = 1.0, 
         legend.lab = NULL,
        legend.line = 2 
           ){
  
   info<- legendLayout
   
   image.plot( legend.only=TRUE,
              add=TRUE,
              smallplot = info$smallplot,
              col= col,
              zlim = zlim,
              legend.lab = legend.lab,
              legend.line= legend.line, 
               axis.args = axis.args,
              legend.args = legend.args,
              legend.cex = legend.cex,
           legend.shrink = info$legend.shrink,
              legend.mar = info$legend.mar,
            legend.width = info$legend.width,
              horizontal = info$horizontal
              )
}
