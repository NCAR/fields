


setupColorBar<- function( horizontal = FALSE,
                        legend.shrink = 0.9, 
                        legend.width = 1.2, 
                        legend.mar = ifelse(horizontal, 3.1, 5.1)
){
  # save current graphics settings
   
  
  
  temp <- imageplot.setup(add = FALSE, 
                          legend.shrink = legend.shrink, 
                          legend.width = legend.width ,
                          legend.mar= legend.mar, 
                          horizontal = horizontal,
                          bigplot = NULL,
                          smallplot = NULL)
  Info<- list( smallplot = temp$smallplot,
               oldPar =  par(no.readonly = TRUE) )
  par(plt = temp$bigplot)
  assign( ".legendInfo",Info, pos=1)
}

addColorBar<- function( legend.lab = NULL,
           legend.line= 2, 
           axis.args = NULL,
           legend.args = NULL,
           legend.cex=1.0, 
           col=NULL,
           zlim=NULL
           ){
  oldPar <- par(no.readonly = TRUE) 
  if( is.null(col)){
    col<- .colorMapInfo$col
  }
  if( is.null( zlim)){
    zlim<- .colorMapInfo$zlim
  }
  image.plot( legend.only=TRUE,
              add=TRUE,
              smallplot = .legendInfo$smallplot,
              col= col,
              zlim = zlim,
              graphics.reset=FALSE, 
              legend.lab = legend.lab,
              legend.line= legend.line, 
              axis.args = axis.args,
              legend.args = legend.args,
              legend.cex=legend.cex)
  
# reset to current plot window  
  par(oldPar)
  mfg.save <- par()$mfg
  par(mfg = mfg.save, new = FALSE)
  invisible()
}
