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

\name{circulantEmbedding}
\alias{sim.rf}
\alias{circulantEmbedding}
\alias{circulantEmbeddingSetup}
\title{
  Efficiently Simulates a Stationary  1 and 2D Gaussian random fields  
}
\description{
Simulates a stationary Gaussian random field on a regular grid with
unit marginal variance. Makes use of the efficient algorithm based on the FFT know as 
circulant embedding.
}
\usage{
sim.rf(obj)
circulantEmbedding(obj)
circulantEmbeddingSetup( 
     grid, M = NULL, cov.function="stationary.cov", cov.args=NULL, ...)

}
\arguments{
\item{obj}{A list (aka covariance object) that includes information about the covariance
function and the grid for evaluation. Usually this is created by a
setup call to Exp.image.cov, stationary.image.cov, matern.image.cov or
other related covariance functions for \code{sim.rf} (See details below.)
or to \code{circulantEmbeddingSetup} for \code{circulantEmbedding} }

\item{grid}{A list describing the regular grid. \code{length(grid)} is the dimension of the field
(1D 2D etc) and each component are the regular locations in that dimension. }

\item{M}{A vector of dimensions to embed the field. Simulation will be exact if each \code{M[i]} is
larger than \code{2*length(grid)}. The default is to choose a power of 2 larger than this
minima bound.}

\item{cov.function}{A text string with the name of the stationary covariance function to use.
Default is \code{stationary.cov} and general function that takes advantage of some efficiency
in finding distances.}

\item{cov.args}{A list of arguments to include with the covariance function, Eg. theta and 
smoothness for the Matern.}

\item{\dots}{For convenience any other arguments to pass to the covariance function. }
}

\value{
\strong{sim.rf}:
A matrix with the random field values.

\strong{circulantEmbedding}: An array according to the grid values specified in the setup.

\strong{circulantEmbeddingetup}: A list with components 

\code{ "m"    "grid" "dx"   "M"    "wght" "call"}

With the information needed to simulate the field.

}
\details{
The functions \code{circulantEmbedding} and \code{circulantEmbeddingSetup} are more recent \code{fields}
functions, more easy to read, and recommended over \code{sim.rf}. \code{sim.rf} is limited to 2D fields
while \code{circulantEmbedding} can handle any number of dimensions and has some shortcuts to be efficient for the 2D case. 

The simulated field has the marginal variance that is determined by
the covariance function for zero distance. Within fields the
exponential and matern set this equal to one ( e.g. Matern(0) ==1) so
that one simulates a random field with a marginal variance of one. For
stationary.cov the marginal variance is whatever \code{Covariance(0)}  evaluates to and we
recommend that alternative covariance functions also be normalized so
that this is one.

Of course if one requires a Gaussian field with different marginal
variance one can simply scale the result of this function. See the
third example below. 

Both \code{sim.rf} and \code{circulantEmbedding}
take an object that includes some preliminary
calculations and so is more efficient for simulating more than one
field from the same covariance. 

The algorithm using an FFT
known as circulant embedding, may not always work if the correlation
range is large. Specifically the weight function obtained from the FFT of the covariance field will have some negative values. 
A simple fix is to increase the size of the domain 
so that the correlation scale becomes smaller relative to the extent
of the domain. Increasing the size can be computationally expensive,
however, and so this method has some limitations. But when it works it is
an exact simulation of the random field. 

For a stationary model the covariance object ( or list) for \code{circulantEmbedding} should have minmally, the components:
That is 
\code{names( obj)}
should give
\code{ "m" "grid"    "M"    "wght" }

where \code{m} is the number of grid points in each dimension, \code{grid} is a list
with components giving the grid points in each coordinate.  
\code{M} is the size of the larger grid that is used for "embedding" and 
simulation. Usually \code{M = 2*m}  and results in an exact
simulation of the stationary Gaussian field. The default  if \code{M} is not passed
is to find the smallest power of 2 greater than 2*m.  \code{wght} is an array from
the FFT of the covariance function with dimensions \code{M}.
Keep in mind that for the final results only the array
that is within the indices   \code{1: m[i]}  for each dimension \code{i} is retained.
This can give a much larger intermediate array, however, in the computation. 
E.g. if \code{m[1] = 100} and \code{m[2]=200} by default then \code{M[1] = 256} 
and  \code{M[2] = 512}. A 256 X 512 array 
is simluated with to get the 100 by 200 result. 

The easiest way to create the
object for simulation is to use \code{circulantEmbeddingSetup}. 

For the older function \code{sim.rf} one uses the image based covariance functions with \code{setup=TRUE} to create the list for simulation.
See the example below for this usage.


The classic reference for this algorithm is 
Wood, A.T.A. and Chan, G. (1994).
    Simulation of Stationary Gaussian Processes in [0,1]^d . Journal of
Computational and Graphical Statistics, 3, 409-432. Micheal Stein and
Tilman Gneiting have also made some additional contributions to the
algortihms and theory.

}
\seealso{
 \link{stationary.cov},  \link{stationary.image.cov}
}
\examples{

#Simulate a Gaussian random field with an exponential covariance function,  
#range parameter = 2.0 and the domain is  [0,5]X [0,5] evaluating the 
#field at a 100X100 grid.  
  grid<- list( x= seq( 0,5,,100), y= seq(0,5,,100)) 
  obj<- circulantEmbeddingSetup( grid, Covariance="Exponential", theta=.5)
  set.seed( 223)
  look<-  circulantEmbedding( obj)
# Now simulate another ... 
  look2<- circulantEmbedding( obj)
# take a look at first two  
 set.panel(2,1)
 image.plot( grid[[1]], grid[[2]], look) 
 title("simulated gaussian fields")
 image.plot( grid[[1]], grid[[2]], look2) 
 title("another realization ...")
 
# Suppose one requires an exponential, range = 2
# but marginal variance = 10 ( rho in fields notation)
look3<- sqrt( 10)*circulantEmbedding( obj)

\dontrun{
# an interesting 3D field

grid<- list(  1:40,  1:40, 1:16  )

obj<- circulantEmbeddingSetup( grid,
                         cov.args=list( Covariance="Matern", theta=2, smoothness=1.0)
                         )
# NOTE: choice of theta is close to giving a negative weight array
set.seed( 122)
look<- circulantEmbedding( obj )
# look at slices in the 3rd dimension 
set.panel( 4,4)
zr<- range( look)
par( mar=c(1,1,0,0))
for(  k in 1:16){
image( grid[[1]], grid[[2]], look[,,k], zlim= zr, col=tim.colors(256),
       axes=FALSE, xlab="", ylab="")
}


}


# same as first example using the older sim.rf

grid<- list( x= seq( 0,10,length.out=100) , y= seq( 0,10,length.out=100) )
obj<-Exp.image.cov( grid=grid, theta=.75, setup=TRUE)
set.seed( 223)
look<- sim.rf( obj)
# Now simulate another ... 
look2<- sim.rf( obj)

 
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
