%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2018
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

\name{bplot.xy}
\alias{bplot.xy}
\title{
  Boxplots for conditional distribution 
}
\description{
Draws boxplots for y by binning on x. This gives a coarse, but quick,  
representation
of the conditional distrubtion of [Y|X] in terms of boxplots.  
}
\usage{
bplot.xy(x, y, N = 10, breaks = pretty(x, N, eps.correct = 1), plot=TRUE,
         ...)
}
\arguments{
\item{x}{
Vector to use for bin membership
}
\item{y}{
Vector to use for constructing boxplot statistics. 
}
\item{N}{
Number of bins on x. Default is 10.
}
\item{breaks}{
Break points defining bin boundaries. These can be unequally spaced. 
}

\item{plot}{
If FALSE just returns a list with the statistics used for plotting the
box plots, bin centers, etc. --  More stuff than you can imagine!
}
\item{\dots }{
Any other optional arguments passed to the standard \code{boxplot} function. 
}
}
\seealso{
bplot, draw.bplot
}
\examples{
# condition on swim times to see how run times vary
bplot.xy( minitri$swim, minitri$run, N=5)

# bivariate normal corr= .8
set.seed( 123)
x<-rnorm( 2000)
y<- .8*x +  sqrt( 1- .8**2)*rnorm( 200)
#
bplot.xy(x,y)
#
bplot.xy( x,y, breaks=seq( -3, 3,,25) ,
                xlim =c(-4,4), ylim =c(-4,4), col="grey80", lwd=2)
points( x,y,col=3, cex=.5)
}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
