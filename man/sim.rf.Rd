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

\name{sim.rf}
\alias{sim.rf}
\title{
  Simulates a Stationary Gaussian random field  
}
\description{
Simulates a stationary Gaussian random field on a regular grid with
unit marginal variance.
}
\usage{
sim.rf(obj)
}
\arguments{

\item{obj}{
A covariance object that includes information about the covariance
function and the grid for evaluation. Usually this is created by a
setup call to Exp.image.cov, stationary.image.cov, matern.image.cov or
other related covariance functions. (See details below.)
}

\item{\dots}{
Additional arguments passed to a particular method.}

}
\value{
A matrix with the random field values 
}
\details{
The simulated field has the marginal variance that is determined by
the covariance function for zero distance. Within fields the
exponential and matern set this equal to one ( e.g. Matern(0) ==1) so
that one simulates a random field with a marginal variance of one. For
stationary.cov the marginal variance is \code{cov.function(0)} and we
recommend that alternative covariance functions also be normalized so
that this is one.

Of course if one requires a Gaussian field with different marginal
variance one can simply scale the result of this function. See the
third example below. 

This function takes an object that includes some preliminary
calculations and so is more efficient for simulating more than one
field from the same covariance. However, the algorithm using a 2-d FFT
(known as circulant embedding) may not always work if the correlation
range is large.  The simple fix is to increase the size of the domain
so that the correlation scale becomes smaller relative to the extent
of the domain. Increasing the size can be computationally expensive
however and so this method has some limitations. But when it works it is
and exact simulation of the random field. 

For a stationary model the covariance object should have  the components: 

names( obj) 
 "m"    "n"    "grid" "N"    "M"    "wght",   

where m and n are the number of grid points in x and y, grid is a list
with components x and y giving the grid points in each coordinate.  N
and M is the size of the larger grid that is used for
simulation. Usually M = 2*m and N =2*n and results in an exact
simulation of the stationary Gaussian field.  wght is a matrix from
the FFT of the covariance function.  The easiest way to create this
object is to use for example Exp.image.cov with setup=T ( see below).

The classic reference for this algorithm is 
Wood, A.T.A. and Chan, G. (1994).
    Simulation of Stationary Gaussian Processes in [0,1]^d . Journal of
Computational and Graphical Statistics, 3, 409-432. Micheal Stein and
Tilman Gneiting have also made some additional contributions to the
algortihms and theory.

}
\seealso{
Exp.image.cov, matern.image.cov, stationary.image.cov 
}
\examples{
#Simulate a Gaussian random field with an exponential covariance function,  
#range parameter = 2.0 and the domain is  [0,5]X [0,5] evaluating the 
#field at a 100X100 grid.  
grid<- list( x= seq( 0,5,,100), y= seq(0,5,,100)) 
obj<-Exp.image.cov( grid=grid, theta=.5, setup=TRUE)
look<- sim.rf( obj)
# Now simulate another ... 
look2<- sim.rf( obj)

# Suppose one requires an exponential, range = 2
# but marginal variance = 10 ( rho in fields notation)
look3<- sqrt( 10)*  sim.rf( obj)

# take a look at first two  
set.panel(2,1)
 image.plot( grid$x, grid$y, look) 
 title("simulated gaussian fields")
 image.plot( grid$x, grid$y, look2) 
 title("another realization ...")
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
