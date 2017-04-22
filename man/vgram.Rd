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

\name{vgram}
\alias{vgram}
\alias{crossCoVGram}
\alias{plot.vgram}
\alias{boxplotVGram}
\alias{getVGMean}
\title{
  Traditional or robust variogram methods for spatial data
}
\description{
\code{vgram} computes pairwise squared differences as a function of distance. 
Returns an S3 object of class "vgram" with either raw values or statistics from 
binning. \code{crossCoVGram} is the same as \code{vgram} but differences are 
taken across different variables rather than the same variable.

\code{plot.vgram} and \code{boxplotVGram} create lineplots and boxplots of 
vgram objects output by the \code{vgram} function.  \code{boxplotVGram} plots 
the base R boxplot function, and plots estimates of the mean over the boxplot.

The \code{getVGMean} function returns the bin centers and means of the \code{vgram} 
object based on the bin breaks provided by the user.
}
\usage{
vgram(loc, y, id = NULL, d = NULL, lon.lat = FALSE, 
                    dmax = NULL, N = NULL, breaks = NULL, 
                    type=c("variogram", "covariogram", "correlogram"))

crossCoVGram(loc1, loc2, y1, y2, id = NULL, d = NULL, lon.lat = FALSE, 
                        dmax = NULL, N = NULL, breaks = NULL, 
                        type=c("cross-covariogram", "cross-correlogram"))

boxplotVGram(x, N=10, breaks = pretty(x$d, N, eps.correct = 1), plot=TRUE, plot.args, ...)

\method{plot}{vgram}(x, N=10, breaks = pretty(x$d, N, eps.correct = 1), add=FALSE, ...)

getVGMean(x, N = 10, breaks = pretty(x$d, N, eps.correct = 1)) 
}
\arguments{
\item{loc}{
Matrix where each row is the coordinates of an observed point
of the field
}
\item{y}{
Value of the field at locations
}
\item{loc1}{
Matrix where each row is the coordinates of an observed point
of field 1
}
\item{loc2}{
Matrix where each row is the coordinates of an observed point
of field 2
}
\item{y1}{
Value of field 1 at locations
}
\item{y2}{
Value of field 2 at locations
}
\item{id}{
A 2 column matrix that specifies which variogram differnces to find. 
If omitted all possible pairing are found. 
This can used if the data has an additional covariate that determines
proximity, for example a time window. 
}
\item{d}{
Distances among pairs indexed by id. If not included distances from from
directly from loc. 
}
\item{lon.lat }{
If true, locations are assumed to be longitudes and latitudes
and distances found are great circle distances (in miles see 
\link{rdist.earth}).  Default is FALSE. 
}
\item{dmax}{
Maximum distance to compute variogram.
}
\item{N}{
Number of bins to use. The break points are found by the \code{pretty} function and so ther may not be exactly N bins. Specify the breaks explicity if you want excalty N bins. 
}
\item{breaks}{
Bin boundaries for binning variogram values. Need not be equally spaced
but must be ordered. 
}
\item{x}{
An object of class "vgram" (an object returned by \code{vgram})
}
\item{add}{
If \code{TRUE}, adds empirical variogram lineplot to current plot.  Otherwise 
creates new plot with empirical variogram lineplot.
}
\item{plot}{
If \code{TRUE}, creates a plot, otherwise returns variogram statistics output by 
\code{bplot.xy}.
}
\item{plot.args}{
Additional arguments to be passed to \code{plot.vgram}.
}
\item{type}{
One of "variogram", "covariogram", "correlogram", "cross-covariogram", and 
"cross-correlogram".  \code{vgram} supports the first three of these and 
\code{crossCoVGram} supports the last two.
}
\item{...}{
Additional argument passed to \code{plot} for \code{plot.vgram} or to 
\code{bplot.xy} for \code{boxplotVGram}.
}
}
\value{
\code{vgram} and \code{crossCoVGram} return a "vgram" object containing the 
following values:
\item{vgram}{Variogram or covariogram values} 
\item{d}{Pairwise distances} 
\item{call}{Calling string} 
\item{stats}{Matrix of statistics for values in each bin. 
  Rows are the summaries returned by the stats function or describe.  
  If not either breaks or N arguments are not supplied then this
  component is not computed.} 
\item{centers}{Bin centers.}

If \code{boxplotVGram} is called with \code{plot=FALSE}, it returns a 
list with the same components as returned by \code{bplot.xy}
}

\section{References}{
See any standard reference on spatial statistics. For example 
Cressie, Spatial Statistics  
}
\author{John Paige, Doug Nychka}
\seealso{
\link{vgram.matrix}, \link{bplot.xy}, \link{bplot}
}
\examples{
#
# compute variogram for the midwest ozone field  day 16
# (BTW this looks a bit strange!)
#
data( ozone2)
good<- !is.na(ozone2$y[16,])
x<- ozone2$lon.lat[good,] 
y<- ozone2$y[16,good]

look<-vgram( x,y, N=15, lon.lat=TRUE) # locations are in lon/lat so use right
#distance
# take a look:
plot(look, pch=19)
#lines(look$centers, look$stats["mean",], col=4)

brk<- seq( 0, 250,, (25 + 1) ) # will give 25 bins.
 
## or some boxplot bin summaries

boxplotVGram(look, breaks=brk, plot.args=list(type="o"))
plot(look, add=TRUE, breaks=brk, col=4)

#
# compute equivalent covariogram, but leave out the boxplots
#
look<-vgram( x,y, N=15, lon.lat=TRUE, type="covariogram")
plot(look, breaks=brk, col=4)

#
# compute equivalent cross-covariogram of the data with itself 
#(it should look almost exactly the same as the covariogram of 
#the original data, except with a few more points in the 
#smallest distance boxplot and points are double counted)
#
look = crossCoVGram(x, x, y, y, N=15, lon.lat=TRUE, type="cross-covariogram")
plot(look, breaks=brk, col=4)

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
