bubblePlot= function( x, y, z, col=tim.colors(256),
                      horizontal = FALSE,
                      legend.shrink = 0.9, 
                      legend.width = 1.2, 
                      legend.mar = ifelse(horizontal, 3.1, 5.1),
                      size=1.0,add=FALSE,
                      ...){
  # NOTE: need to use
  # when just adding points to a plot 
  #  add=TRUE
  # use setupLegend before plotting to 
  # include a legend color scales
  
  # save current graphics settings
  old.par <- par(no.readonly = TRUE)
  # adjust values of x,y, and z if (x,y) passed in first argument
   x<- as.matrix( x)
  if( dim(x)[2]==2){
    z<- y 
    y<- x[,2]
    x<- x[,1]
  }
# color table for z  values  
  ctab= color.scale( z, col)
  if(!add){
    setupLegend( 
      legend.shrink = legend.shrink, 
      legend.width = legend.width ,
      legend.mar = legend.mar, 
      horizontal = horizontal
       )  
      plot( x,y, type="n", ...) 
  } 
# this is it  
  points( x,y, col=ctab, pch=16, cex=size)
  
# add legend if setup has been called   
  if(exists(".legendInfo")){
    addLegend(col=col, zlim=range(z) )
  }
  
}

