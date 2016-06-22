%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2016
%# University Corporation for Atmospheric Research (UCAR)
%# Contact: Douglas Nychka, nychka@ucar.edu,
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

\name{Chicago ozone test data}
\alias{ChicagoO3}
\alias{ozone}
\title{
  Data set of ozone measurements at 20 Chicago monitoring stations. 
}
\description{
This data set used be named  \code{ozone} but was 
changed to avoid conflict with other packages.
The \code{ChicagoO3} data is a list of 
components, x and y.  x 
component is longitude and latitude position of each of the 20 
Chicago monitoring stations, y is the average 
daily ozone values over the time period 6/3/87-8/30/87.
These data are used extensively for the test scripts and 
simple examples. The lasting scientific value is probably 
minimal. 

}

\format{
This data set is a list containing the following components: 

\describe{
 \item{lon.lat}{ Longitude-latitude positions of monitoring stations.  }
 \item{x}{An approximate Cartesian set of coordinates for the locations 
where the units are in miles. The origin is in the center of the 
locations. }
 \item{y}{ Average daily ozone values over 1987 summer. }
 }
}
\source{
AIRS, the EPA air quality data base. 
}
\seealso{
Tps, Krig  
}
\examples{
fit<- Tps(ChicagoO3$x, ChicagoO3$y) 
# fitting a surface to ozone measurements. 
surface( fit, type="I")
}
\keyword{datasets}
% docclass is function
% Converted by Sd2Rd version 1.21.
