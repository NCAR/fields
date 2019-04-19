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
\name{Krig.replicates}
\alias{Krig.replicates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Collapse repeated spatial locations into unique locations
}
\description{
In case that several observations are available for a single spatial location find the 
group means and replicate variability
}
\usage{
Krig.replicates(out, x, y, Z, weights=rep( 1, length(y)), verbose = FALSE)
}

\arguments{
  \item{out}{ A list with components  \code{x}, \code{y}, \code{weights},  and  possibily 
\code{Z}.}
\item{x}{Spatial locations.}
\item{y}{Spatial observations}
\item{Z}{Spatial covariates.}
\item{weights}{Weights proportional to reciprocal varainces of observations.}

\item{verbose}{ If TRUE print out details for debugging. }
}
\details{
This function figures out which locations are the same and within the function fast.1way 
use \code{tapply} to find replicate group means and standard deviations. 
NOTE: it is assumed the Z covariates are unique at the locations. Currently these functions
can not handle a model with common spatial locations but different values for the Z covariates.

}
\value{
A list with components:
  \item{yM }{Data at unique locations and where more than one observation is
 available this is the mean of the replicates.}

  \item{xM }{Unique spatial locations.}
\item{weightsM}{Weights matching the unique lcoations proportional to reciprocal variances 
This is found as a combination of the original weights at each location.}
\item{ZM}{Values of the covariates at the unique lcoations.}
\item{uniquerows}{Index for unique rows of \code{x}.}
\item{shat.rep, shat.pure.error}{Standard deviation of pure error estimate based on replicate groups
 (and adjusting for possibly different weights.)}
\item{rep.info}{Integer tags indicating replicate groups.}

}

\author{
Douglas Nychka
}

\examples{

#create  some spatial replicates
 set.seed( 123)
 x0<- matrix( runif(10*2), 10,2)
 x<-  x0[ c(rep(1,3), 2:8, rep( 9,5),10) , ]
 y<-  rnorm( 16)
 
 out<- Krig.replicates( x=x, y=y)
# compare 
# out$yM[1] ;  mean( y[1:3])
# out$yM[9] ; mean( y[11:15])
# mean( y[ out$rep.info==9])
 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
