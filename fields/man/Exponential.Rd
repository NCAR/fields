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

\name{Exponential, Matern, Radial Basis}
\alias{Exponential}
\alias{Matern}
\alias{Matern.cor.to.range}
\alias{RadialBasis}
\title{Covariance functions}
\description{
Functional form of covariance function assuming the argument is a
distance between locations.  As they are defined here, they are in 
fact correlation functions.  To set the marginal variance (sill) 
parameter, use the \code{sigma} argument in \code{mKrig} or \code{Krig}.  
To set the nugget variance, use te \code{tau2} argument in 
\code{mKrig} or \code{Krig}.
}
\usage{
Exponential(d, range = 1, alpha = 1/range, phi=1.0)
Matern(d , range = 1,alpha=1/range, smoothness = 0.5, 
       nu= smoothness, phi=1.0) 
Matern.cor.to.range(d, nu, cor.target=.5, guess=NULL,...)
RadialBasis(d,M,dimension, derivative = 0)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{d}{ Vector of distances or for \code{Matern.cor.to.range} just a single distance. }

  \item{range}{ Range parameter default is one. Note 
       that the scale can also be specified  through the "aRange" 
           scaling argument  used in fields covariance functions) }

  \item{alpha}{1/range }
  
  \item{phi}{This parameter option is added to be compatible with older
  versions of fields and refers to the marginal variance of the process.
  e.g. \code{phi* exp( -d/aRange)} is the exponential covariance for points
  separated by distance and range aRange. Throughout fields this parameter
  is equivalent to sigma and it recommended that sigma be used. If one is
  simulating random fields. See the help on \code{\link{sim.rf}} for
  more details. } 

  \item{smoothness}{ Smoothness parameter in Matern. Controls the number
of derivatives in the process. Default is 1/2 corresponding to an exponential 
covariance.}

  \item{nu}{ Same as smoothness}
  \item{M}{Interpreted as a spline M is the order of the derivatives in the
    penalty.}
  \item{dimension}{Dimension of function} 
   \item{cor.target}{Correlation used to match the range parameter. Default is .5.}
   \item{guess}{An optional starting guess for solution. This should not be needed.}
   \item{derivative}{If greater than zero finds the first derivative of this function.}
   \item{\dots}{Additional arguments to pass to the bisection search function.}
}

\details{
Exponential: 

exp( -d/range)

Matern:

   con*(d\^nu) * besselK(d , nu )

 Matern covariance function transcribed from Stein's book page 31
nu==smoothness, alpha ==  1/range

 GeoR parameters map to kappa==smoothness and phi == range
check for negative distances

\code{con} is a constant that normalizes the expression to be 1.0 when d=0. 

Matern.cor.to.range: 
    This function is useful to find  Matern covariance parameters that are 
comparable for different smoothness parameters. Given a distance \code{d}, 
smoothness \code{nu},  target correlation \code{cor.target} and
range \code{aRange}, this function determines numerically the value of 
aRange so that

\code{Matern( d, range=aRange, nu=nu) == cor.target}

See the example for how this might be used.

Radial basis functions:
\preformatted{
   C.m,d  r**(2m-d)        d- odd

   C.m,d  r**(2m-d)ln(r)    d-even
}
where C.m.d is a constant based on spline theory and r is the radial distance
between points. See \code{radbas.constant} for the computation of the constant.
NOTE: Earlier versions of fields used ln(r^2) instead of ln(r) and so differ by a factor of 2. 
}
\value{

For the covariance functions: a vector of covariances.

For Matern.cor.to.range: the value of the range parameter.  

}
\references{ Stein, M.L. (1999) Statistical Interpolation of Spatial Data: Some Theory for Kriging. Springer, New York.}
\author{Doug Nychka}
\seealso{stationary.cov, stationary.image.cov, Wendland,stationary.taper.cov
  rad.cov}

\examples{
# a Matern correlation function 
 d<- seq( 0,10,,200)
 y<- Matern( d, range=1.5, smoothness=1.0)
 plot( d,y, type="l")

# Several Materns of different smoothness with a similar correlation 
# range

# find ranges for nu = .5, 1.0 and 2.0 
# where the correlation drops to .1 at a distance of 10 units.

 r1<- Matern.cor.to.range( 10, nu=.5, cor.target=.1)
 r2<- Matern.cor.to.range( 10, nu=1.0, cor.target=.1)
 r3<- Matern.cor.to.range( 10, nu=2.0, cor.target=.1)

# note that these equivalent ranges
# with respect to this correlation length are quite different
# due the different smoothness parameters. 

 d<- seq( 0, 15,,200)
 y<- cbind(  Matern( d, range=r1, nu=.5),
             Matern( d, range=r2, nu=1.0),
             Matern( d, range=r3, nu=2.0))

 matplot( d, y, type="l", lty=1, lwd=2)
 xline( 10)
 yline( .1)
}
\keyword{spatial}% at least one, from doc/KEYWORDS
