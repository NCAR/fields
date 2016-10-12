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
MLESpatialProcess(x, y, weights = rep(1, nrow(x)), Z = NULL, mKrig.args
                 = NULL, cov.function = "stationary.cov", cov.args =
                 list(Covariance = "Matern", smoothness = 1),
                 lambda.guess = 0.5, theta.start = NULL, theta.range =
                 NULL, gridN = 20, optim.args = NULL, na.rm = TRUE,
                 verbose = FALSE, ...)
                           
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
A matrix of spatial locations with rows indexing location
  and columns the dimension (e.g. longitude/latitude)
}
\item{y}{
Spatial observations
}
\item{lambda.start}{
The initial guess for lambda, the nugget to sill ratio.
}
\item{theta.start}{
The initial guess for theta, the correlation range parameter.
}
\item{cov.function}{
The name of the covariance function (See help on Krig for details. )
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
\item{Distance}{
Distance function to use in covariance. Used for setting initial guess 
  and range for theta.
}
\item{verbose}{
If TRUE print out intermediate information for debugging. 
}
\item{\dots}{
Additional arguments to pass to the mKrig function.
}
}
\details{
	MLESpatialProcess is designed to be a simple and easy to use function to 
	maximize the likelihood for a Gaussian spatial process. For certain fixed,
	covariance parameters, the likelihood is maximized over the nugget and sill
	parameters using the \code{mKrig} function.  \code{lambda} and \code{theta} 
	are optimized using the \code{mKrigMLEJoint} function on a log scale.
	
	MLESpatialProcess.fast uses the \code{optim}
	function to maximize the likelihood computed from the \code{mKrig} function. It is more
	efficient in the computation as it does not find a full eigen decomposition with each new value
	of theta and maximizes the likelihood over theta and lambda simultaneously. 
	
	Note the likelihood can be maximized analytically over the  parameters of the fixed part and with the nugget (sigma) and sill (rho)  reduced to the single parameter lambda= sigma^2/rho. So fixing any other covariance parameters the likelihood is maximzed numerically over lambda and theta. The differences between these two functions is due to the differences between the definition of the
	restricted likelihood used in Krig and the conventional likelihood used in mKrig.
	 
  
}
\value{
	\code{MLESpatialProcess}:
	A list that includes components:
\code{theta.MLE, rho.MLE, sigma.MLE, lambda.MLE} being the maximum
 likelihood estimates of these
parameters. The component \code{REML.grid} is a two column matrix
 with the
first column being the theta grid and the second column being the
 profiled and restricted  likelihood for that value of theta. Here profile means that
  the likelihood has already been evaluated at the maximum over sigma
   and rho for this value of theta.
	      \code{eval.grid} is a more complete "capture" of the 
	      evaluations being a 
	      6 column matrix with the parameters theta, lambda, sigma,
	       rho, profile likelihood and the effective degrees of
	        freedom. This is just last row of 
	      \code{lambda.est} returned by the core function \code{Krig}
	      
	     \code{MLESpatialProcess.fast} here the returned value is limited because this
	     function isbuilt around calls to \code{mKrig}. Returned value is a list with components:
	     \code{pars}, the MLEs for theta, rho, sigma and lambda,  
		  \code{logLikelihood},values of the log likelihood at the maximum, 
		     \code{eval.grid}, a table with the results from evaluating different combinations of 
		     parameters,

            \code{converge}, convergence flag from optim (0=Successfull) and number of evaluations used to find maximum.
            and \code{call}, the calling arguments. 
            
            }

\author{
Doug Nychka, John Paige
}

\seealso{
\code{\link{Krig}}, \code{\link{mKrigMLEGrid}}, \code{\link{mKrigMLEJoint}}, \code{\link{optim}}, \code{\link{fastTps.MLE}}, \code{\link{spatialProcess}}
}
\examples{
#
#

#generate observation locations
n=200
x = matrix(runif(2*n), nrow=n)

#generate observations at the locations

trueTheta = .2
trueLambda = .1
Sigma = exp( -rdist(x,x) /trueTheta ) 
# y = t(chol(Sigma))%*% (rnorm(n))  +  trueLambda* rnorm( n)
y = t(chol(Sigma))\%*\% (rnorm(n))  +  trueLambda* rnorm( n)

#Use exponential covariance, assume the true range parameter is known
out = MLESpatialProcess(x, y, 
                        cov.args=list(Covariance="Exponential", range=trueTheta)
                        )

#Use exponential covariance, use a range to determine MLE of range parameter
\dontrun{

#Use Matern covariance, compute joint MLE of range, smoothness, and lambda.  
#This may take a few seconds
testSmoothness = c(.5, 1, 2)
for( nu in testSmoothness){
  out = MLESpatialProcess(x, y, cov.args=list(Covariance="Matern"), smoothness=nu) 
  print( out$summary)
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
obj2$summary
obj3$summary

  }
 \dontrun{
# fits for first 10 days from ozone data
data( ozone2) 	 
NDays<- 10
O3MLE<- matrix( NA, nrow= NDays, ncol=7)
for( day in 1: NDays){
	cat( day, " ")
	ind<- !is.na(ozone2$y[day,] )
	x<- ozone2$lon.lat[ind,]
	y<- ozone2$y[day,ind]
	print( length( y))
	O3MLE[day,]<- MLESpatialProcess( x,y,
	            Distance="rdist.earth")$MLEJoint$summary
}
# NOTE: names of summary:
#[1] "lnProfileLike.FULL" "lambda"            
#[3] "theta"              "sigmaMLE"          
#[5] "rhoMLE"             "funEval"           
#[7] "gradEval" 
plot( log(O3MLE[,2]), log(O3MLE[,3]))
}  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
