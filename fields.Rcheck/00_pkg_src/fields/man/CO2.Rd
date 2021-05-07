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
\name{CO2}
\alias{CO2}
\alias{CO2.true}
\docType{data}
\title{Simulated global CO2 observations}
\description{
 This is an example of moderately large spatial data set and consists of 
simulated CO2 concentrations that are irregularly sampled from a lon/lat 
grid. Also included is the complete CO2 field  (CO2.true) used to generate the 
synthetic observations.}
 \usage{data(CO2)  }
\format{
  The format of \code{CO2} is a list with two components:
\itemize{
\item lon.lat: 26633x2 matrix of the longitude/latitude locations. 
 These are a subset of  a larger lon/lat grid (see example below).
\item y:   26633 CO2 concentrations in parts per million.
}
 The format of \code{CO2.true}
 is a list in  "image" format with components:
\itemize{
 \item x longitude grid values.
 \item y latitude grid values.
 \item z an image matrix with CO2 concentration in parts per million
 \item mask a logical image that indicates with grid locations were 
selected for the synthetic data set \code{CO2}.
}

}
\details{ 
This data was generously provided by Dorit Hammerling and Randy Kawa as a 
test example for the spatial analysis of remotely sensed (i.e. 
satellite) and irregular observations. The synthetic data is based on a 
true CO2 field simulated from a geophysical, numerical model. 
}
\examples{
\dontrun{

data(CO2)
#
# A quick look at the observations with world map
quilt.plot( CO2$lon.lat, CO2$y)
world( add=TRUE)

# Note high concentrations in Borneo (biomass burning), Amazonia and
# ... Michigan (???).

# spatial smoothing using the wendland compactly supported covariance
# see help( fastTps) for details
# First smooth using locations and Euclidean distances 
# note taper is in units of degrees 
out<-fastTps( CO2$lon.lat, CO2$y, aRange=4, lambda=2.0) 
#summary of fit note about 7300 degrees of freedom 
# associated with fitted surface
 print( out)
# image plot on a grid  (this takes a while)
surface( out, type="I", nx=300, ny=150)
# smooth with respect to great circle distance 
out2<-fastTps( CO2$lon.lat, CO2$y, lon.lat=TRUE,lambda=1.5, aRange=4*68) 
print(out2)
#surface( out2, type="I", nx=300, ny=150)

# these data are actually subsampled from a grid. 
# create the image object that holds the data
#

temp<- matrix( NA, ncol=ncol(CO2.true$z), nrow=nrow(CO2.true$z))
temp[ CO2.true$mask] <- CO2$y

# look at gridded object. 
 image.plot(CO2.true$x,CO2.true$y, temp)

# to predict _exactly_ on this grid for the second fit;
# (this take a while)
look<- predictSurface( out2, grid.list=list( x=CO2.true$x, y=CO2.true$y))
image.plot(look)

}
}
\keyword{datasets}
