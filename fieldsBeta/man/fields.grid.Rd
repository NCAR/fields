%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2018
%# University Corporation for Atmospheric Research (UCAR)
%# Contact: Douglas Nychka, nychka@mines.edu,
%# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
%#
%# This program is free software; you can redistribute it and/or modify
%# it under the terms of the GNU General Public License as published by
%# the Free Software Foundation; either version 2 of the License, or
%# (at your option) any later version.
%# This program is distributed in the hope that it will be useful,
%# but WITHOUT ANY WARRANTY; without even the implied warranty of
%# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%# GNU General Public License for more details.
%#
%# You should have received a copy of the GNU General Public License
%# along with the R software environment if not, write to the Free Software
%# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%# or see http://www.r-project.org/Licenses/GPL-2    

\name{fields.grid}
\alias{mKrig.grid}
\title{
Using MKrig for predicting on a grid.
}

\description{
This is an extended example for using the sparse/fast interpolation
methods in mKrig to evaluate a Kriging estimate on a large grid. 
}

\details{
\code{mKrig} is a flexible function for surface fitting  using 
a spatial process model. It can also exploit sparse matrix methods forlarge data sets by using a compactly supported covariance. 
The example below shows how ot evaluate a solution on a big grid. (Thanks to Jan Klennin for  this example.)
}
\examples{
x<- RMprecip$x
y<- RMprecip$y

Tps( x,y)-> obj

# make up an 80X80 grid that has ranges of observations
# use same coordinate names as the x matrix

glist<- fields.x.to.grid(x, nx=80, ny=80) # this is a cute way to get a default grid that covers x

# convert grid list to actual x and y values ( try plot( Bigx, pch="."))
    make.surface.grid(glist)-> Bigx 

# include actual x locations along with grid. 
    Bigx<- rbind( x, Bigx)

# evaluate the surface on this set of points (exactly)

    predict(obj, x= Bigx)-> Bigy

# set the range for the compact covariance function 
# this will involve  less than 20 nearest neighbors that have
# nonzero covariance
# 

     V<- diag(c( 2.5*(glist$lon[2]-glist$lon[1]), 
                 2.5*(glist$lat[2]-glist$lat[1])))
\dontrun{
# this is an interplotation of the values using a compact 
# but thin plate spline like covariance. 
    mKrig( Bigx,Bigy, cov.function="wendland.cov",k=4, V=V, 
                 lambda=0)->out2 
# the big evaluation this takes about 45 seconds on a Mac G4 latop
    predictSurface( out2, nx=400, ny=400)-> look
}

# the nice surface
\dontrun{    
    surface( look)
    US( add=TRUE, col="white")
}

}
 

\keyword{hplot}
