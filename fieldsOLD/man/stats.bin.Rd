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

\name{stats.bin}
\alias{stats.bin}
\title{
  Bins data and finds some summary statistics. 
}
\description{
Cuts up a numeric vector based on binning by a covariate and applies the
fields 
stats function to each group 
}
\usage{
stats.bin(x, y, N = 10, breaks = NULL)
}
\arguments{
\item{x}{
Values to use to decide bin membership
}
\item{y}{
A vector of data
}
\item{N}{
Number of bins. If the breaks is missing there are N bins equally spaced
on the range of x.
}
\item{breaks}{
The bin boundaries. If there are N+1 of these there will be N bins. 
The bin widths can be unequal.
}
}
\value{
A list with several components. stats is a matrix with columns indexing
the bins and
rows being summary statistics found by the stats function. These are:
 number of obs, mean, sd, min, quartiles, max and number of NA's. 
(If there is no data for a given bin, NA's are filled in. )
breaks are the breaks passed to the function and centers are the bin
centers. 
}
\seealso{
bplot, stats
}
\examples{
u<- rnorm( 2000)
v<- rnorm( 2000)
x<- u
y<- .7*u + sqrt(1-.7**2)*v

look<- stats.bin( x,y) 
look$stats["Std.Dev.",]

data( ozone2)
# make up a variogram day 16 of Midwest daily ozone ...
look<- vgram( ozone2$lon.lat, c(ozone2$y[16,]), lon.lat=TRUE)

# break points
brk<- seq( 0, 250,,40)

out<-stats.bin( look$d, look$vgram, breaks=brk)
# plot bin means, and some quantiles  Q1, median, Q3
matplot( out$centers, t(out$stats[ c("mean", "median","Q1", "Q3"),]), 
type="l",lty=c(1,2,2,2), col=c(3,4,3,4), ylab="ozone PPB")
}
\keyword{univar}
% docclass is function
% Converted by Sd2Rd version 1.21.
