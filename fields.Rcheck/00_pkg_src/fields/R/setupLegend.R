


setupLegend<- function( horizontal = FALSE,
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

