colorBar<- function(breaks,
                    smallplot,
                    colorTable,
                    horizontal = FALSE,
                    lab.breaks,
                    axis.args,
                    legend.lab,
                   legend.line = 2,
                    legend.args,
                    legend.cex = 1,
                    lowerTriangle = FALSE,
                    upperTriangle = NULL
                    ){
  # Following code draws the legend using the image function
  # and a one column image.
  # What might be confusing is the values of the "image" are the same 
  # as the locations on the legend axis.
  # Moreover the image values are in the middle of each breakpoint category
  # thanks to Tobias Nanu Frechen and Matthew Flickinger 
  # for sorting out some problems with the breaks position in the legend.
  ix <- 1:2
  # if triangles for end points  upper or lower limits in breaks should 
  # be omitted and the colorTable also amended. 
  lowerColor<- NULL
  upperColor<- NULL
  if( lowerTriangle){
    breaks<- breaks[-1]
    lowerColor<- colorTable[1]
    colorTable<- colorTable[-1]
  }
  if( upperTriangle){
    n1<- length(breaks)
    breaks<- breaks[-n1]
    n2<- length(colorTable)
    upperColor<- colorTable[n2]
    colorTable<- colorTable[-n2]
  }
  iy<- breaks
  N<- length(colorTable)
  nBreaks<- length( breaks)
  midpoints<- (breaks[1:(nBreaks-1)] +  breaks[2:nBreaks] )/2
  iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints)) 
  #
  # next par call sets up a new plotting region just for the legend strip
  # at the smallplot coordinates
  par(new = TRUE, pty = "m", plt = smallplot, err = -1)
  # draw color scales the two  cases are horizontal/vertical 
  # add a label if this is passed.
  if (!horizontal) {
    image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
          ylab = "", col = colorTable, breaks=breaks)
  }
  else {
    image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
          ylab = "", col = colorTable, breaks=breaks)
  }
 # add triangles at ends  
  addTriangle<- !( is.null(lowerColor)& is.null(upperColor))
  if(addTriangle){
    addColorBarTriangle(
                         lowerColor=lowerColor,
                         upperColor=upperColor,
                         horizontal=horizontal
    )
    }
  
  # create the argument list to draw the axis
  # this avoids 4 separate calls to axis and allows passing extra
  # arguments.
  if (!is.null(lab.breaks)) {
    # axis with labels at break points
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
                        at = breaks, labels = lab.breaks), axis.args)
  }
  else {
    # If lab.breaks is not specified ( with or without breaks), pretty
    # tick mark locations and labels are computed internally,
    # or as specified in axis.args at the function call
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
                   axis.args)
  }
  #
  # now add the axis to the legend strip.
  # notice how all the information is in the list axis.args
  do.call("axis", axis.args)
  #
  # add a label to the axis if information has been  supplied
  # using the mtext function. The arguments to mtext are
  # passed as a list like the drill for axis (see above)
  #
  if (!is.null(legend.lab)) {
    legend.args <- list(text = legend.lab,
                        side = ifelse(horizontal, 1, 4), 
                        line = legend.line, # this may need to be tuned
                         cex = legend.cex)
  }
  # add the label using mtext function
  if (!is.null(legend.args)) {
    do.call(mtext, legend.args)
  }
  #
}


