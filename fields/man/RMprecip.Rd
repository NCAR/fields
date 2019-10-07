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

\name{RMprecip}
\alias{RMprecip}
\alias{RMelevation}
\alias{PRISMelevation}
\title{
Monthly total precipitation (mm) for August 1997 in the Rocky Mountain 
Region and some gridded 4km elevation data sets (m).  
}
\description{
\code{RMprecip} is a useful spatial data set of moderate size consisting of 806
locations. 
\code{PRISMelevation} and \code{RMelevation} are gridded elevations for the 
continental US and Rocky Mountain region at 4km resolution. 
Note that the gridded elevations from the PRISM data product are
different than the exact station elevations. (See example below.)
}

\format{
The data set  \code{RMprecip} is a list containing the following components: 

\describe{
\item{x}{
Longitude-latitude position of monitoring stations. Rows names are station id codes consistent
with the US Cooperative observer network.
The ranges for these coordinates are  [-111, -99] for longitude and [35,45] for latitude.
}
\item{elev}{
Station elevation in meters. 
}
\item{y}{
Monthly total precipitation in millimeters.
for August, 1997
}
}

The data sets 
 \code{PRISMelevation} and  
 \code{RMelevation} are lists 
in the usual R grid format for images and contouring

They have  the following components:
\describe{
\item{x}{
Longitude grid at approximately 4km resolution}
\item{y}{
Latitude grid at approximately 4km resolution}
\item{z}{
Average  elevation for grid cell in meters
}

}

These elevations and the companion grid formed the basis for the
103-Year High-Resolution Precipitation Climate Data Set for the
Conterminous United States
( see \url{http://www.prism.oregonstate.edu/documents/PRISM_downloads_FTP.pdf}
and also archived at the National
Climate Data Center. This work is authored by Chris Daly
\url{http://www.prism.oregonstate.edu} and his PRISM group but had some
contribution from the Geophysical Statistics Project at NCAR  
and is an interpolation of the observational data to a 4km grid that
takes into account topography such as elevation and aspect.

}
\details{
Contact Doug Nychka for the binary file \code{RData.USmonthlyMet.bin}
and  information on its source. 

\preformatted{
# explicit source code to create the RMprecip data
dir <- "" # include path to data file 
load(paste(dir, "RData.USmonthlyMet.bin", sep="/")
#year.id<-  1963- 1895
year.id<- 103
#pptAUG63<- USppt[ year.id,8,]
loc<- cbind(USpinfo$lon, USpinfo$lat)
xr<- c(-111, -99)
yr<- c( 35, 45)
station.subset<-  (loc[,1]>= xr[1]) & (loc[,1] <= xr[2]) & (loc[,2]>= yr[1]) & (loc[,2]<= yr[2])
ydata<-  USppt[ year.id,8,station.subset]
ydata <- ydata*10 #  cm -> mm conversion
xdata<- loc[station.subset,]
dimnames(xdata)<- list( USpinfo$station.id[station.subset], c( "lon", "lat"))
xdata<- data.frame( xdata)
good<- !is.na(ydata)
ydata<- ydata[good]
xdata<- xdata[good,]
     
test.for.zero.flag<- 1
test.for.zero( unlist(RMprecip$x), unlist(xdata), tag="locations")
test.for.zero( ydata, RMprecip$y, "values")
}
}
\examples{
# this data set was created  the 
# historical data  taken from 
# Observed monthly precipitation, min and max temperatures for the coterminous US 
# 1895-1997
# NCAR_pinfill 
# see the Geophysical Statistics Project datasets page for the supporting functions 
# and details. 

# plot 
quilt.plot(RMprecip$x, RMprecip$y)
US( add=TRUE, col=2, lty=2)

# comparison of station elevations with PRISM gridded values

data(RMelevation)

interp.surface( RMelevation, RMprecip$x)-> test.elev

plot( RMprecip$elev, test.elev, xlab="Station elevation", 
ylab="Interpolation from PRISM grid")
abline( 0,1,col="blue")

# some differences  with high elevations probably due to complex
# topography!

#
# view of Rockies looking from theSoutheast

save.par<- par(no.readonly=TRUE)

par( mar=c(0,0,0,0))

# fancy use of persp with shading and lighting.
persp( RMelevation, theta=75, phi= 15, 
          box=FALSE, axes=FALSE, xlab="", ylab="", 
         border=NA,
         shade=.95, lphi= 10, ltheta=80,
         col= "wheat4", 
         scale=FALSE, expand=.00025)

# reset graphics parameters and a more conventional image plot.
par( save.par)
image.plot(RMelevation, col=topo.colors(256))
US( add=TRUE, col="grey", lwd=2)
title("PRISM elevations (m)")
}
\keyword{datasets}
% docclass is data
% Converted by Sd2Rd version 1.21.
