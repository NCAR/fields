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
\name{NorthAmericanRainfall}
\alias{NorthAmericanRainfall}
\docType{data}
\title{
Observed North American summer precipitation from the
historical climate network.
}
\description{
Average rainfall in tenths of millimeters for the months
of June, July and August for the period 1950-2010. Data is based on 
1720 stations located in North America.
}
%\usage{data(NorthAmericanRainfall)}
\format{
  The format is a list with components:
"longitude" "latitude"  "precip"    "elevation" "precipSE" "trend" "trendSE" "type" "x.s" "sProjection"
with elevation in meters, longitude as (-180,180), latitude as (-90, 90) and precipitaion in 1/10 mm 
( precip/254  converts to inches of rainfall)

\code{precip} is the intercept for 1980.5 when a straight line least squares regression is fit to 
each station's record. SE is the companion standard error from the least squares fit.
If the station is complete, then \code{precip} and \code{precipSE} will just be the mean and standard deviation adjusted for a linear trend. The estimated trend  \code{trend} and and its standard error \code{trendSE} are also included. 
 Also due to the centering, for complete data the intercept and trend estimate will be uncorrelated. The component \code{type} indicates whether the station has been "adjusted" (see below) or is still in
  "unadjusted" form. 

\code{x.s} is a useful transformation of locations into stereographic coordinates that reduces the
inflation of North Canada due to the usual lon/lat coordinates. Specifically it is found by: 
\preformatted{
  library(mapproj)
  xStereo<- mapproject( NorthAmericanRainfall$lon,NorthAmericanRainfall$lat, projection="stereographic")
  NorthAmericanRainfall$x.s<- cbind( xStereo$x, xStereo$y)
  NorthAmericanRainfall$projection<- .Last.projection
}
Use \code{NorthAmericanRainfall$orientation} to access the stereographic projection orientation.

}

\source{
The monthly data used to construct this summary was generously provided by Xuebin Zhang, however,
the orignal source is freely available as the Global Historical Climate Network Version 2 Precipitation
quality controlled, curated and served by the US National Center for Environmental Information. 
The adjusted data from this archive has been modified from its raw form to make the record more homogenous. Heterogenities can come from a variety of sources such as a moving the station a short distance or changes in instruments. See \url{https://www.ncdc.noaa.gov} and goto data-access and then GHCN. 
}
\examples{
data(NorthAmericanRainfall)
x<- cbind(NorthAmericanRainfall$longitude,  NorthAmericanRainfall$latitude)
y<- NorthAmericanRainfall$precip
quilt.plot( x,y)
world( add=TRUE)

Zstat<- NorthAmericanRainfall$trend / NorthAmericanRainfall$trendSE
quilt.plot( x, Zstat)

}
\keyword{datasets}
