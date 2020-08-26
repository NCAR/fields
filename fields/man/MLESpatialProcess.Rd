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
\name{MLESpatialProcess}
\alias{MLESpatialProcess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Estimates key covariance parameters for a spatial process.
%%  ~~function to do ... ~~
}
\description{
Maximizes the likelihood to determine the nugget variance (sigma^2), the sill
( rho) and the range (theta) for a spatial process.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MLESpatialProcess( x, y, weights = rep(1, nrow(x)), Z = NULL,
                            mKrig.args = NULL,
                          cov.function = "stationary.cov", 
                              cov.args = list(Covariance = "Matern",
                                            smoothness = 1),  
			     gridTheta = NULL,
                                 gridN = 15,
                            optim.args = NULL,
                                 na.rm = TRUE,
                               verbose = FALSE,
                               abstol  = 1e-4,
                                  REML = FALSE,
                      cov.params.start = NULL,
                               ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
A matrix of spatial locations with rows indexing location
  and columns the dimension (e.g. longitude/latitude)
}
\item{y}{
The spatial observations. This can be a matrix of reoplicated spatial data. 
}
  \item{weights}{
Precision ( 1/variance) of each observation
}
  \item{Z}{
Linear covariates to be included in fixed part of the 
  model that are distinct from the default low order 
  polynomial in \code{x}
}
\item{mKrig.args}{A list containing other objects to pass to mKrig.}


\item{gridN}{Number of points to use in grid search over theta.}
\item{cov.function}{
The name of the covariance function (See help on Krig for details. )
}

\item{gridTheta}{ A grid of range parameters to search over. The default is to use a range based on the
quantiles of pairwise distances of the spatial locations.
}

\item{cov.args}{
A list with arguments for the covariance functions. These are 
  usually parameters and other options such as the type of 
  distance function. 
}


\item{optim.args}{
Additional arguments passed to the optim function for likelihood 
  maximization.  The default value is:
  \code{optim.args = list(method = "BFGS", 
                    control = list(fnscale = -1, parscale = c(0.5, 0.5), 
                                   ndeps = c(0.05,0.05)))}
}
\item{na.rm}{If TRUE remove missing values in y and corresponding locations in x.}

\item{verbose}{If TRUE print out intermediate information for debugging.}

\item{abstol}{Absolute tolerance used to judeg convergence in optim.}

\item{REML}{If TRUE use maximize the restricted Likelihood instead of the concentrated
likelihood.(Preliminary experience suggests this does not make much difference.) }

\item{cov.params.start}{A list with each component being the name of a covariance paramater to
optimize over and the value being the starting value. This is the same format used by optim. If
NULL no additional parameters are optimized and their values are fixed in the cov.args list.}

\item{\dots}{
Additional arguments to pass to the mKrig function.
}
}
\details{
	MLESpatialProcess is designed to be a simple and easy to use function for
	maximizing the likelihood for a Gaussian spatial process. For other  fixed,
	covariance parameters, the likelihood is maximized over the nugget and sill
	parameters using the \code{mKrig} function.  \code{lambda} and \code{theta} 
	are optimized using the \code{mKrigMLEJoint} function on a log scale.  This is
	a wrapper for two lower level functions to make \code{spatialProcess} readable. We
	recommend using   \code{mKrigMLEJoint} and \code{mKrigMLEGrid} for direct
	optimzation and grid searches over parameters and consulting
	those two commented functions for an exact description of how
	the computations are being done. At the very basic
	level the likelihood function is evaluated by calling the
	\code{mKrig} function with the covariance parameters being
	adjusted to test different values. So the lowest level
	likelihood evaluate and what one might use at the highest
	level are the same workhorse function, \code{mKrig}. 
	
	Note the likelihood can be maximized analytically over the 
	parameters of the fixed
	part of the spatial model and with the nugget (sigma) and
	sill (rho)  reduced to the
	single parameter lambda= sigma^2/rho.  The
likelihood is maximized numerically over lambda and theta if
there are additional
covariance parameters ( such as smoothness for the Matern)
these need to be fixed
and so the MLE is found for the covariance conditional on these
additional
parameter values. From a practical point of view it is often
difficult to estimate
just these three from a moderate spatial data set and the user
is encourage to try
different combinations of fixing covariance parameters with ML
for the remaining
ones.  
	
MLESpatialProcess.fast is an older fields function also using
the \code{optim}
function to maximize the likelihood computed from the 
\code{mKrig} function. It will
eventually be removed from later versions of fields but perhaps
is
still useful as a cross
check on these newer functions  
	 
}
\value{
%list( summary= MLEJoint$summary, MLEGrid= MLEGrid, MLEJoint=MLEJoint, 
%          MLEProfileLambda=MLEProfileLambda, call=match.call(), 
%         timingTable= timingTable)
	\code{MLESpatialProcess}:
A list (with sublists!) that includes components:
	
\code{summary} A vector with the optimized parameters (aka MLEs) along with log likelihood 
at the maximum and information from optim on convergence. Also included are the effecive 
degrees of freedom of the model at the MLEs and the GCV function value. 

\code{MLEJoint}  A list with components: 
\code{summary, lnLike.eval,  optimResults, pars.MLE,    parTransform}.
\code{summary} gives the results from the \strong{joint}
optimization over 
all parameters. \code{lnLike.eval} is a table 
that has captured all the likelihood evaluations as optim has
searched for the maxmimum. Note that
if the method has converged many of these parameter value
choices will be close the final converged
value but often still useful to investigate the likelihood
surface around the maximum.  \code{optimResults} is the list
returned by the optim function. 

\code{pars.MLE} are the parameters specifically optimized over.
\code{par.Transform} the function used to
transform the parameters for a more suitable optimization.
The default is a  log transformation for
lambda and theta. 

\code{MLEGrid} A list with results from profiling the likelihood over theta, the range parameter. 
Here  the component \code{summary} is a table giving the values tried for theta and the corresponding parameters that
maximize the likelihood with this fixed value for theta. Also included in the table are some convergence
information, effective degrees of freedom, and the GCV criterion. 

\code{MLEProfileLambda} Results from profiling the likelihood over lambda. Format is the same as that described for \code{MLEGrid}.

\code{call} The call made to this function.

\code{timingTable} Some timing information for the joint optimzation, as the column labeled timeOptim, and
the two profiling computations with timeGrid being for theta and timeProfile being for lambda. 

	     \code{MLESpatialProcess.fast} has been depreciated and is included for backward compatibility.  
            
            }

\author{
Doug Nychka, John Paige
}

\seealso{
\code{\link{Krig}}, \code{\link{mKrigMLEGrid}}, \code{\link{mKrigMLEJoint}}, \code{\link{optim}}, \code{\link{fastTps.MLE}}, \code{\link{spatialProcess}}
}
\examples{
#
# Trying out on a simulated data set
# Generate observation locations (50 is really too small but this is
#  just to make this run quickly)

  n <- 50
 
  set.seed(124)
  x = matrix(runif(2*n), nrow=n)
#generate observations at the locations
  trueTheta = .1
  trueSigma = .01
  covMat = exp( -rdist(x,x) /trueTheta )
# As rcode:    y = t(covMat)%*% (rnorm(n))  +  trueSigma * rnorm( n)
  y  <-  t(covMat) \%*\% (rnorm(n))  +  trueSigma * rnorm( n)
# Use exponential covariance estimate constant function for mean
  out <-  MLESpatialProcess(x, y, 
                        smoothness = .5,
                        mKrig.args = list( m = 1), 
                        )
  out$summary
\dontrun{
#Now a (small) grid search over the smoothness
# add replicated data  fields ( i.e. independent copies drawn from same covariance model)
# to make this stable
  set.seed(124)
# Y = t(covMat)%*% matrix(rnorm(n*200), n,200)  +  trueSigma * matrix((rnorm(n*200)), n,200)
  Y = t(covMat)\%*\% matrix(rnorm(n*200), n,200)  +
                  trueSigma * matrix((rnorm(n*200)), n,200)
#This may take a few seconds
testSmoothness = c(.5, 1.5, 2.0)
for( nu in testSmoothness){
  out = MLESpatialProcess(x, Y, cov.args = list(Covariance="Matern"),
                  smoothness = nu, cov.params.start= list( lambda=.5)) 
  print( out$MLEJoint$summary)
}

}

# example with a covariate  
\dontrun{
data(COmonthlyMet)
ind<- !is.na( CO.tmean.MAM.climate)
x<- CO.loc[ind,]
y<- CO.tmean.MAM.climate[ind]
elev<- CO.elev[ind]
obj2<-  MLESpatialProcess( x,y)
obj3<-  MLESpatialProcess( x,y, Z=elev)

# elevation makes a difference
obj2$MLEJoint$summary
obj3$MLEJoint$summary

  }
 \dontrun{
# fits for first 10 days from ozone data
data( ozone2) 	 
NDays<- 5
O3MLE<- NULL
for( day in 1: NDays){
	cat( day, " ")
	ind<- !is.na(ozone2$y[day,] )
	x<- ozone2$lon.lat[ind,]
	y<- ozone2$y[day,ind]
	print( length( y))
	 out<- MLESpatialProcess( x,y,
	            Distance="rdist.earth")$MLEJoint$summary
	O3MLE<- rbind( O3MLE, out)
}
# NOTE: names from summary:
#[1] "lnProfileLike.FULL" "lambda"             "theta"             
#[4] "sigmaMLE"           "rhoMLE"             "funEval"           
#[7] "gradEval"           "eff.df"             "GCV"  
plot( log(O3MLE[,"lambda"]), log(O3MLE[,"theta"]))
}  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
