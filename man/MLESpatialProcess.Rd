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
MLESpatialProcess(x, y, theta.grid=NULL, par.grid=NULL, lambda.grid=NULL, 
                  cov.function = "stationary.cov", 
                  cov.args = list(Covariance = "Matern", smoothness = 1), 
                  optim.args = NULL, ngrid = 10, niter = 15, tol = 0.01, 
                  Distance = "rdist", nstep.cv = 50, verbose = FALSE, 
                  doMKrig=FALSE, ...)
                           
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
\item{theta.grid}{
Grid of theta parameter values to use for grid search in maximizing 
  the Likelilood. The defualt is do an initial grid search on ngrid 
  points with the range at the 3 an d 97 quantiles of the pairwise 
  distances.  If only two points are passed then this is used as the 
  range for a sequence of ngrid points.  Note that this is only used 
  when \code{doMKrig==FALSE}.
}
\item{par.grid}{
Grid list of covariance parameters and the values to use in the grid 
  search for maximizing the Likelilood. All combinations of parameter 
  values are used in the grid search.  If only two values 
  for a parameter are passed then this is used as the range for a 
  sequence of ngrid points.  If a full set of parameter values 
  is passed, it is recommended they be distributed on a log scale.  
  Make sure to put all fixed parameters in cov.args rather than 
  par.grid.  Note that this is only used when \code{doMKrig==TRUE}.
}
\item{lambda.grid}{
Grid list of lambda values to use for grid search in maximizing 
  the Likelilood. If NULL, automatically set to \code{10^seq(-6, 1, by=1)}  
  If only two points are passed then this is used as the range 
  for a sequence of ngrid points.  If a full grid.list of 
  parameter values is passed, it is recommended they be distributed 
  on a log scale.  If lambda is fixed, put the fixed value as a ... 
  argument rather than in lambda.grid.  Note that this is only used when 
  \code{doMKrig==TRUE}
}
\item{cov.function}{
The name of the covariance function (See help on Krig for details. )
}
\item{cov.args}{
A list with arguments for the covariance functions. These are 
  usually parameters and other options such as the type of 
  distance function. 
}
\item{ngrid}{
Number of points in grid search over parameters. 
}
\item{nstep.cv}{
Number of grid points to use in GCV or REML coarse
	search for optimum.  Note that this is only used 
  when \code{doMKrig==FALSE}.
}
\item{optim.args}{
Additional arguments that would also be included in calls to the optim 
  function in the final joint likelihood maximization with initial 
  lambda and covariance guesses set to the Tps maximum.  The default 
  value in this function is:
  \code{optim.args = list(method = "BFGS", 
                        control=list(fnscale = -1, 
                                     ndeps = rep(log(1.1), 2), 
                                     reltol=1e-02, maxit=3))}
  Note that this argument is only used when \code{doMKrig==TRUE}.
}
\item{niter}{
Max number of iterations for the golden section search to maximize over theta.
    Note that this is only used when \code{doMKrig==FALSE}.
}
\item{tol}{
Tolerance to declare convergence. 
}
\item{doMKrig}{
If \code{TRUE}, uses \code{mKrig}.  If \code{FALSE}, uses \code{Krig}.  This 
  will also change which input arguments are used.
}
\item{Distance}{
Distance function to use in covariance.
}
\item{verbose}{
If TRUE print out intermediate information for debugging. 
}
\item{\dots}{
Additional arguments to pass to the Krig or mKrig function depending on 
  which is used.
}
}
\details{
	MLESpatialProcess is designed to be a robust but perhaps slow function to 
	maximize the likelihood for a Gaussian spatial process. For certain fixed,
	covariance parameters, the likelihood is maximized over the nugget and sill
	parameters using the \code{Krig} or \code{mKrig} function. An outer 
  optimization finds the maximum over other specified covariance parameters 
  as well. See the help(Krig) for details of the restricted maximum 
  likelihood criterion (REML).
	
	MLESpatialProcess.fast uses the \code{optim}
	function to maximize the likelihood computed from the \code{mKrig} function. It is more
	efficient in the computation as it does not find a full eigen decomposition with each new value
	of theta and maximizes the likelihood over theta and lambda simultaneously. 
	
	Note the likelihood can be maximized analytically over the  parameters of the fixed part and with the nugget (sigma) and sill (rho)  reduced to the single parameter lambda= sigma^2/rho. So fixing any other covariance parameters the likelihood is maximzed numerically over lambda and theta. The differences between these two functions is due to the differences between the definition of the
	restricted likelihood used in Krig and the conventional likelihood used in mKrig.
	 
  In general, it is recommended to perform joint optimization using 
  \code{mKrig}, which evaluates the log-likelihood over the grid of 
  lambda and covariance parameter values, interpolating them with a 
  thin-plate spline.  Afterwards, a final joint optimization is 
  performed using \code{mKrig.MLE.joint} with the initial guess set 
  to the thin-plate spline maximum.  It may be more advantageous to 
  use \code{Krig}, however, when only performing optimization over 
  lambda.
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
\code{\link{Krig}}, \code{\link{mKrig.MLE}}, \code{\link{mKrig.MLE.joint}}, \code{\link{optim}}, \code{\link{fastTps.MLE}}, \code{\link{spatialProcess}}
}
\examples{
#
# examples with doMKrig==TRUE
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
                        cov.args=list(Covariance="Exponential", range=trueTheta), 
                        doMKrig=TRUE)

#Use exponential covariance, use a range to determine MLE of range parameter
\dontrun{
testThetas = seq(from=trueTheta/2, to=2*trueTheta, length=6)
par.grid=list(theta=testThetas)
out = MLESpatialProcess(x, y, 
                        cov.args=list(Covariance="Exponential"), 
                        par.grid=par.grid, doMKrig=TRUE)                        

#Use exponential covariance, use a range to determine MLE of range 
#parameter, set custom lambda.grid
testLambdas= seq(from=trueLambda/2, to=2*trueLambda, length=6)
out = MLESpatialProcess(x, y, cov.args=list(Covariance="Exponential"), 
                        lambda.grid=testLambdas, par.grid=par.grid, doMKrig=TRUE)

#Use Matern covariance, compute joint MLE of range, smoothness, and lambda.  
#This may take a few seconds
testSmoothness = c(.5, 1, 2)
par.grid=list(range=testThetas, smoothness=testSmoothness)
out = MLESpatialProcess(x, y, cov.args=list(Covariance="Matern"), 
                        par.grid=par.grid, doMKrig=TRUE)
}

#
# examples with doMKrig==FALSE
#
	N<- 100
  set.seed(123)
  x<- matrix(runif(2*N), N,2)
  theta<- .2
  Sigma<-  Matern( rdist(x,x)/theta , smoothness=1.0)
  Sigma.5<- chol( Sigma)
  sigma<- .1
#  F.true<- t( Sigma.5)%*%  rnorm(N)
  F.true<- t( Sigma.5) \%*\%  rnorm(N)
  Y<-  F.true +  sigma*rnorm(N)
# find MLE for  sigma rho and theta  smoothness fixed  first
# data set
  obj<- MLESpatialProcess( x,Y)
  obj$pars
  # profile likelihood over theta
  plot(obj$eval.grid[,1], obj$eval.grid[,6], xlab="theta",
  ylab= "log Profile likelihood", type="p" )
  xline( obj$pars["theta"], col="red")
# log likelihood surface over theta and  log lambda
image.plot( obj$logLikelihoodSurface$x,
               obj$logLikelihoodSurface$y, obj$logLikelihoodSurface$z,
               xlab="theta (range)", ylab="log lambda" )
# MLE               
  points( obj$pars[1], log(obj$pars[2]), pch=16, col="magenta", cex=1.2)      

# using "fast" version  
obj.fast<- MLESpatialProcess.fast( x,Y)
obj.fast$pars 
# points where likelihood evaluated:
quilt.plot( log( obj.fast$eval.grid[,1:2] ), obj.fast$eval.grid[,7],
                   xlab="log(theta)",ylab="log(lambda)")

# parameters are slightly different due to the differences of REML and the full likelihood.
                 
# example with a covariate  
\dontrun{
data(COmonthlyMet)
obj2<-  MLESpatialProcess( CO.loc, CO.tmean.MAM.climate)
obj3<-  MLESpatialProcess( CO.loc, CO.tmean.MAM.climate, Z= CO.elev)
ind<- !is.na( CO.tmean.MAM.climate)
obj4<-  MLESpatialProcess.fast( CO.loc[ind,], CO.tmean.MAM.climate[ind],
               Z= CO.elev[ind])
# elevation makes a difference
obj2$pars
obj3$pars
obj4$pars
  }
 \dontrun{
# fits for ozone data
data( ozone2) 	 
NDays<- nrow( ozone2$y)
O3MLE<- matrix( NA, nrow= NDays, ncol=4)
dimnames(O3MLE)<- list(NULL, c("theta",  "lambda", "rho",    "sigma"))
for( day in 1: NDays){
	cat( day, " ")
	O3MLE[day,]<- MLESpatialProcess( ozone2$lon.lat, ozone2$y[day,],
	 Distance="rdist.earth")$pars
}
plot( log(O3MLE[,1]), log(O3MLE[,3]))
}  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
