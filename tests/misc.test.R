 # fields is a package for analysis of spatial data written for
  # the R software environment .
  # Copyright (C) 2018
  # University Corporation for Atmospheric Research (UCAR)
  # Contact: Douglas Nychka, nychka@ucar.edu,
  # National Center for Atmospheric Research,
  # PO Box 3000, Boulder, CO 80307-3000
  #
  # This program is free software; you can redistribute it and/or modify
  # it under the terms of the GNU General Public License as published by
  # the Free Software Foundation; either version 2 of the License, or
  # (at your option) any later version.
  # This program is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  # GNU General Public License for more details.

suppressMessages(library(fields))
options( echo=FALSE)
set.seed( 234)
test.for.zero.flag<-1

y<- runif( 30)
lev<- sort(sample( 1:5,30, replace=TRUE))
w<- runif( 30)*.1+1
y<- as.matrix(y)

# compute by loop
hold<- rep( NA, 5)
for( k in 1:5){
  ind<- lev==k
  hold[k]<- sum( y[ind,]*w[ind])/ sum( w[ind])}

look<- fast.1way( lev, y, w)
test.for.zero( look$means, hold, tag="fast.1way means")

# now vectorized case

ytemp<- cbind( y, y-10, y+10)
look2<- fast.1way( lev, ytemp, w)
test.for.zero( look2$means[,2], hold-10, tag="fast.1way vectorized means")


cat("All done with testing misc functions", fill=TRUE)
options(echo=TRUE)
