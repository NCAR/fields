envelopePlot <- function(x1, y1, x2 = x1, y2,
                         col ="thistle1" , lineCol = "thistle3", ...) {
#  sort the curves -- just in case they are passed out of order
    ind<- order( x1)
    x1<- x1[ind]
    y1<- y1[ind]
    ind<- order( x2)
    x2<- x2[ind]
    y2<- y2[ind]
   
  polygon(c(x1, rev(x2)), c(y1, rev(y2)), col = col, border = NA, ...)
  lines(x1, y1, lwd = 3, col = lineCol)
  lines(x2, y2, lwd = 3, col = lineCol)
}
