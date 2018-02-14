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
\name{REML.test}
\Rdversion{1.1}
\alias{REML.test}
\alias{MLE.Matern}
\alias{MLE.Matern.fast}
\alias{MLE.objective.fn}
\alias{MaternGLS.test}
\alias{MaternGLSProfile.test}
\alias{MaternQR.test}
\alias{MaternQRProfile.test}
\title{
Maximum Likelihood estimates for some Matern covariance parameters. 
}
\description{
For a fixed smoothness (shape) parameter these functions provide
different ways of estimating and testing restricted and profile
likehiloods for the Martern covariance parameters. \code{MLE.Matern}
is a simple function that finds the restricted maximum likelihood
(REML) estimates of the sill, nugget and range parameters (\code{rho,
sigma2 and theta}) of the Matern covariance functions.  The remaining
functions are primarily for testing.
}
\usage{


MLE.Matern(x, y, smoothness, theta.grid = NULL, ngrid = 20,
                 verbose = FALSE, niter = 25, tol = 1e-05,
                 Distance = "rdist", m = 2, Dmax = NULL, ...)

MLE.Matern.fast(x, y, smoothness, theta.grid = NULL, ngrid=20, verbose=FALSE,
                                         m=2, ...)
MLE.objective.fn( ltheta,info, value=TRUE)

MaternGLSProfile.test(x, y, smoothness = 1.5, init = log(c(0.05,1)))
MaternGLS.test(x, y, smoothness = 1.5, init = log(c(1, 0.2, 0.1)))
MaternQR.test (x, y, smoothness = 1.5, init = log(c(1, 0.2, 0.1))) 
MaternQRProfile.test (x, y, smoothness = 1.5, init = log(c(1))) 

REML.test(x, y, rho, sigma2, theta, nu = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{Dmax}{ Maximum distance for grid used to evaluate the fitted
  covariance function.}

  \item{Distance}{Distance function used in finding covariance.}

  \item{x}{ A matrix of spatial locations with rows indexing location
    and columns the dimension (e.g. longitude/latitude)}
  
\item{y}{ Spatial observations}

  \item{smoothness}{Value of the Matern shape parameter.}

  \item{theta.grid}{ Grid of theta parameter values to use for grid
  search in maximizing the Likelilood. The defualt is do an initial
  grid search on ngrid points with the range at the 3 an d 97
  quantiles of the pairwise distances.If only two points are passed
  then this is used as the range for a sequence of ngrid points.}

   \item{ngrid}{Number of points in grid search.}

\item{init}{Initial values of the parameters for optimization. For the
first three functions these are in the order rho, theta sigma2 and in
a log scale. For MaternQRProfile.test initial value is just
log(theta). }

\item{verbose}{If TRUE prints more information.}

\item{rho}{ Marginal variance of Matern process (the "sill") }

\item{sigma2}{Variance of measurement error (the "nugget")}

\item{theta}{Scale parameter (the "range")}

\item{nu}{Smoothness parameter}


\item{ltheta}{ log of range parameter}

\item{info}{A list with components \code{x,y, smoothness, ngrid} that
pass the information to the optimizer. See details below.}

\item{value}{If TRUE only reports minus log Profile likelihood with
profile on the range parameter. If FALSE returns a list of
information.}

\item{m}{Polynomial of degree (m-1) will be included in model as a
fixed part.}

\item{niter}{Maximum number of interations in golden section search.}

\item{tol}{Tolerance for convergence in golden section search.}

\item{\dots}{Additional arguments that are passed to the Krig function
in evaluating the profile likelihood.}

}

\details{ 
 \code{MLE.Matern} is a simple function to find the maximum likelihood
estimates of using the restricted and profiled likeilihood that is
intrinsic to the ccomputations in \code{Krig}. The idea is that the
likelihood is concentrated to the parameters lambda and theta. (where
lambda = sigma2/rho). For fixed theta then this is maximized over
lambda using \code{Krig} and thus concetrates the likelihood on
theta. The final maximization over theta is implemented as a golden
section search and assumes a convex function. All that is needed is
for three theta grid points where the middle point has a larger
likelihood than the endpoints. In practice the theta grid defualts to
a 20 points equally spaced between the .03 and .97 quantiles of the
distribution of the pairwise distances. The likelihood is evaluated at
these points and a possible triple is identified.  If no exists from
the grid search the function returns with NAs for the parameter
estimates.  Note that due to the setup of the golden section search
the computation actually minimizes minus the log likelihood.
\code{MLE.Matern.fast} is a similar function but replaces the
optimaiztion step computed by Krig to a tighter set of code in the
function \code{MLE.objective.fn}.  See also \code{mKrigMLEGrid} for an
alternative and streamlined function using \code{mKrig} rather than
\code{Krig}.
}
\value{
For MLE.Matern (and MLE.Matern.fast)

\item{smoothness}{Value of the smoothness function}

\item{pars}{MLE for rho, theta and sigma}

\item{REML}{Value of minus the log restricted Profile likelihood at
the maxmimum}

\item{trA}{Effective degrees of freedom in the predicted surface based
on the MLE parameters.}

\item{REML.grid}{Matrix with values of theta and the log likelihood
from the initial grid search.}

}
\author{
Doug Nychka
}
\note{
See the script REMLest.test.R and Likelihood.test.R in the tests directory to 
see how these functions are used to check the likelihood expressions. 
}
\examples{
# Just look at one day from the ozone2 
data(ozone2)

out<- MLE.Matern( ozone2$lon.lat,ozone2$y[16,],1.5, ngrid=8)
plot( out$REML.grid)
points( out$pars[2], out$REML, cex=2)
xline( out$pars[2], col="blue", lwd=2)
\dontrun{
# to get a finer grid on initial search:
out<- MLE.Matern( ozone2$lon.lat,ozone2$y[16,],1.5,
                      theta.grid=c(.3,2), ngrid=40) 

# simulated data  200 points uniformly distributed
set.seed( 123)
x<- matrix( runif( 2*200), ncol=2)
n<- nrow(x)
rho= 2.0
sigma= .05
theta=.5

Cov.mat<-  rho* Matern( rdist(x,x), smoothness=1.0, range=theta)
A<- chol( Cov.mat)
gtrue<- t(A) \%*\% rnorm(n)
gtrue<- c( gtrue)
err<-  rnorm(n)*sigma
y<- gtrue + err
out0<- MLE.Matern( x,y,smoothness=1.0) # the bullet proof version
# the MLEs and -log likelihood at maximum
print( out0$pars)
print( out0$REML)

out<- MLE.Matern.fast( x,y, smoothness=1.0) # for the impatient
# the MLEs:
print( out$pars) 
print( out$REML)


# MLE for fixed theta (actually the MLE from out0) 
# that uses MLE.objective.fn directly
info<- list( x=x,y=y,smoothness=1.0, ngrid=80)
# the MLEs:
out2<- MLE.objective.fn(log(out0$pars[2]), info, value=FALSE)
print( out2$pars)
}

\dontrun{
# Now back to Midwest ozone pollution ...
# Find the MLEs for ozone data and evaluate the Kriging surface.
  data(ozone2)
  out<- MLE.Matern.fast( ozone2$lon.lat,ozone2$y[16,],1.5)
#use these parameters to fit surface ....
  lambda.MLE<- out$pars[3]/out$pars[1]
  out2<- Krig( ozone2$lon.lat,ozone2$y[16,] , Covariance="Matern",
              theta=out$pars[2], smoothness=1.5, lambda= lambda.MLE)
  surface( out2) # uses default lambda -- which is the right one.

# here is another way to do this where the new lambda is given in 
# the predict step
  out2<- Krig( ozone2$lon.lat,ozone2$y[16,] , Covariance="Matern",
               theta=out$pars[2], smoothness=1.5)
# The default lambda is that found by GCV
# predict on a grid but use the MLE value for lambda:
  out.p<- predictSurface(out2, lambda= lambda.MLE)
  surface(out.p) # same surface!
}

# One could also use mKrig with a fixed lambda to compute the surface. 

\dontrun{
# looping  through all the days of the ozone data set.
  data( ozone2)
  x<- ozone2$lon.lat
  y<- ozone2$y
  out.pars<- matrix( NA, ncol=3, nrow=89)

  for ( k in 1:89){
    hold<- MLE.Matern.fast( x,c(y[k,]), 1.5)$pars
    cat( "day", k," :", hold, fill=TRUE)
    out.pars[k,]<- hold }
}

}
\keyword{spatial}

