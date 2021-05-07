addLegend<- function(
                col = NULL,
               zlim = NULL,
          axis.args = NULL,
        legend.args = NULL,
         legend.cex = 1.0, 
         legend.lab = NULL,
        legend.line = 2 
           ){
  
  if( !exists( ".legendInfo")) {
        stop("setupLegend has not been called")
   }
   legendInfo<-  get(".legendInfo")
   if( is.null(col)){
     if( !exists(  ".colorMap")){
        stop("the argument col  has not been specifed and
                color.scale has not been called to specify one.")
     }  
     colorMapInfo<- get(".colorMapInfo")
     col<- colorMapInfo$col
     zlim<-colorMapInfo$zlim
    }
   image.plot( legend.only=TRUE,
              add=TRUE,
              smallplot = legendInfo$smallplot,
              col= col,
              zlim = zlim,
              legend.lab = legend.lab,
              legend.line= legend.line, 
                axis.args = axis.args,
              legend.args = legend.args,
              legend.cex = legend.cex
              )
}
