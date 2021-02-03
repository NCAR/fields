bubblePlot= function( x, y, z, 
                      col = tim.colors(256),
                      horizontal = FALSE,
                      legend.cex = 1.0,
                      legend.lab = NULL,
                      legend.line= 2, 
                      legend.shrink = 0.9, 
                      legend.width = 1.2, 
                      legend.mar = ifelse(horizontal, 3.1, 5.1),
                      axis.args = NULL, legend.args = NULL,
                      size=1.0,
                      add=FALSE,
                      ...){
  # NOTE: need to use
  # when just adding points to a plot 
  #  add=TRUE
  # use setupLegend before plotting to 
  # include a legend color scales
  
 
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

      plot( x,y, col=ctab,...) 
      
  } 
  else{
 # just add points to an existing plot
  points( x,y, col=ctab, ...)
  }
  
# save current graphics settings  
  big.par <- par(no.readonly = TRUE)
  mfg.save <- par()$mfg
 # print( big.par$plt)
 # print( big.par$usr)
  
# add legend if setup has been called  
  if(exists(".legendInfo")){
    addLegend( 
      col = col,
      zlim = range(z, na.rm=TRUE),
      axis.args = axis.args, 
      legend.args = legend.args,
              legend.cex = legend.cex,
              legend.lab = legend.lab,
              legend.line = legend.line 
              )
  }
# return to graphics setting of the  scatterplot  
# so more good stuff can be added ...
  par(plt = big.par$plt, xpd = FALSE)
  par(mfg = mfg.save, new = FALSE)
  
}

