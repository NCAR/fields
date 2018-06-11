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
test.for.zero.flag<- 1

data(COmonthlyMet)
y<- CO.tmin.MAM.climate
good<- !is.na( y)
y<-y[good]
x<- CO.loc[good,]
Z<- CO.elev[good]
out<- mKrig( x,y, Z=Z,  cov.function="stationary.cov", Covariance="Matern",
                    theta=4.0, smoothness=1.0, lambda=.1)

out2<- Krig( x,y, Z=Z,  cov.function="stationary.cov", Covariance="Matern",
                    theta=4.0, smoothness=1.0, lambda=.1, GCV=TRUE)

test.for.zero( predict( out), predict(out2), tag="Full prediction")
test.for.zero( predict( out, drop.Z=TRUE), predict(out2, drop.Z=TRUE), tag=" prediction dropping Z")

xnew<- CO.loc[!good,]
Znew<-  CO.elev[!good]
temp1<- predict( out, xnew=xnew, drop.Z=TRUE)
temp2<- predict( out2, x=xnew, drop.Z=TRUE)
test.for.zero( temp1,temp2, tag="new x's dropping Z")

temp1<- predict( out, xnew=xnew, Z=Znew)
temp2<- predict( out2, x=xnew, Z=Znew)
test.for.zero( temp1,temp2, tag="new x's new Z's")

temp1<- predictSurface( out, nx=20, ny=20, drop.Z=TRUE, extrap=TRUE)
temp2<- predictSurface( out2, nx=20, ny=20, drop.Z=TRUE, extrap=TRUE)
test.for.zero( temp1$z,temp2$z, tag="predicting on surface with drop.Z")


cat("all done with mKrig Z tests", fill=TRUE)
options( echo=TRUE)

