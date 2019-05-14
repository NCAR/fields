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
\name{spatialProcess}
\alias{spatialProcess}
\alias{plot.spatialProcess}
\alias{print.spatialProcess}
\alias{print.spatialProcessSummary}
\alias{summary.spatialProcess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Estimates a spatial process model. 
%%  ~~function to do ... ~~
}
\description{
	For a given covariance function estimates the nugget (sigma^2) and process variance (rho)  and the
	range parameter (theta) by restricted maximum likelihood and then computes the 
	spatial model with these estimated parameters.
	Other parameters of the covariance 
	are kept fixed and need to be specified.  
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spatialProcess(x, y, weights = rep(1, nrow(x)), Z = NULL, 
mKrig.args = list(m = 2), 
cov.function = "stationary.cov", cov.args = list(Covariance = "Matern", 
        smoothness = 1), theta = NULL, theta.start = NULL, lambda.start = 0.5,
         theta.range = NULL,
    abstol = 1e-04, na.rm = TRUE, verbose = FALSE, REML = FALSE, ...)
    
\method{summary}{spatialProcess}(object,  ...)
\method{print}{spatialProcess}(x, digits = 4, ...)
\method{print}{spatialProcessSummary}(x, digits = 4, ...)
\method{plot}{spatialProcess}(x, digits = 4, which = 1:4, ...)

}

\arguments{
\item{x}{Observation locations}

\item{y}{Observation values}

\item{weights}{Weights for the error term (nugget) in units of reciprocal variance.}

\item{Z}{A matrix of extra covariates for the fixed part of spatial model.
E.g. elevation for fitting climate data over space. }

\item{mKrig.args}{Arguments passed to the mKrig function.}

\item{cov.function}{A character string giving the name of the  covariance function
  	 for the spatial component.}
  	 
 \item{cov.args}{ A list specifying parameters and other components of the covariance function.}
 
 \item{theta}{If not NULL the range parameter for the covariance is fixed at this value. }
 
 \item{theta.start}{Starting value for MLE fitting of the scale (aka range) parameter. If omitted the starting value is taken from a grid search ove theta.}
 
 \item{lambda.start}{Starting value for MLE fitting of the lambda parameter. Note lambda
  is the ratio of the nugget variance to the process variance. In code variables this is
  \code{sigma^2} divided by \code{rho}. }
  \item{ theta.range }{A range for the ML search to estimate theta. Default is based on quantiles of the location pairwise distances.}

 \item{na.rm}{If TRUE NAs are removed from the data. }

% \item{gridN}{Number of grid points for evaluating profile likelihood over
% theta and
% also over lambda (using MLE for theta).}

\item{object}{A spatialProcess object returned from the spatialProcess function.}
 
 %\item{optim.args}{Arguements to be used  in optim for finding joint
 %MLEs for theta
 %and lambda. }
 
 \item{REML}{ If TRUE the parameters are found by restricted maximum likelihood.}
 
 \item{verbose}{If TRUE print out intermediate information for debugging.}

  \item{\dots}{
  	Any other arguments that will be passed to the \code{mKrig} function and interpreted
  	as additional arguments to the covariance function. E.g. \code{smoothness} for the
  	Matern covariance. 
  	}
  	
  \item{abstol}{The absolute tolerance bound used to judge convergence. This is applied
  to the difference in log likelihood values. }	
  
  \item{digits}{Number of significant digits in printed summary}
  
  \item{which}{The vector 1:4 or any subset of 1:4, giving the plots to draw.
                 See the description ofthese plots below.}

}
\details{
	This function makes many choices for the user in terms of defaults and it is 
	important to be aware of these. 
	The spatial model is 
	
	Y.k= P(x.k) + Z(x.k)\%*\%d2 + g(x.k) + e.k

where ".k" means subscripted by k, Y.k is the dependent variable
observed at location x.k. P is a low degree polynomial (default is a
linear function in the spatial coordinates) and Z is a matrix of covariates (optional) that enter as a linear model the fixed part. g is a mean zero,
Gaussian stochastic process with a marginal variance of rho and a
scale (or range) parameter, theta. The measurement errors, e.k, are
assumed to be uncorrelated, normally distributed with mean zero and
standard deviation sigma. If weights are supplied then the variance of e is assumed to be  \code{sigma^2/ weights}.


Perhaps the most important aspect of this  function is that
 	the range (theta), nugget (sigma**2) and process variance (rho) parameters
 	for the covariance are estimated by restricted maximum
 	likelihood and this is the model that is then  used for spatial
 	prediction. Geostatistics usaually refers to sigma**2 + rho as the 
 	"sill" and often these parameters are estimated by variogram fitting rather
 	than maximum likelihood. To be consistent with spline models and to focus
 	on the key part of model we reparametrize as lambda= sigma**2/
 	rho and rho. Thinking about h as the spatial signal and e as the noise lambda
 	can be interpreted  as the noise to signal variance ratio in this spatial
 	context.(See the comparision with fitting the geoR model in the
examples section.)
 	
 	The likelihood and the cross valdiation function
 	can be concentrated to only depend on lambda and theta and so
 	in reported the optimiztation of these two criterion we focus
 	on this form of the parameters.  Once lambda and theta are
 	found, the MLE for rho has a closed form and of course then
 	sigma is then determined from lambda and rho. 
 	
 	Often the lambda
 	parameter is difficult to interpret when covariates and a
 	linear function of the coordinates is included and also when
 	the range becomes large relative to the size of the spatial
 	domain. For this reason it is convenient to report the
 	effective degrees of freedom (also referred to trA in R code and
 	the output summaries) associated with the predicted
 	surface or curve. This measure has a one to one relationship
 	with lamdba and is easier to interpret. For example an eff
 	degrees of freedom that is very small suggests that the
 	surface is rwell represented by a low ordoer
 	polynomial. Degrees of freedom close to the number of
 	locations indicates a surface that is close to interpolating
 	the observations and suggests a small or zero value for the
 	nugget variance.
 	
 	The default covariance model is assumed to follow a Matern
 	 with smoothness set to 1.0.  This is implementd using the
 	 \code{stationary.cov} covariance that can take a argument for
 	 the form of the covariance, a sill and range parameters and
 	 possibily additional parameter might comtrol the shape.
 	 
 	 See the example below how to switch to another model. (Note
 	 that the exponential is also part of the Matern family with
 	 smoothness set to .5. )
  	
 	The parameter estimation is done by \code{MLESpatialProcess}
 	and the returned list from this function is added to the Krig
 	output object that is returned by this function. The estimate
 	is a version of maximum likelihood where the observations are
 	transfromed to remove the fixed linear part of the model.  If
 	the user just wants to fix the range parameter theta then
 	\code{Krig} can be used.
 	
 	NOTE: The defaults for the \code{optim} function used in MLESpatialProcess are:
 	
  \preformatted{
  list(method = "BFGS", 
       control=list(fnscale = -1,
                      ndeps = rep(log(1.1),length(cov.params.start)+1), 
                     abstol = abstol,
                      maxit = 20))
}
                                     

There is always a hazard in providing a simple to use method that
makes many default choices for the spatial model. As in any analysis
be aware of these choices and try alternative models and parameter
values to assess the robustness of your conclusions. Also examine the
residuals to check the adequacy of the fit. See the examples below for
some help in how to do this easily in fields. Also see quilt.plot to
get an quick plot of a spatial field  to discern obvious spatial paterns.


\strong{summary} method forms a list with class \code{spatialProcessSummary} that hasa subset of information from the output object and also creates a table of the estimates of the linear parameters in the fixed part of the model.
With replacated fields there is an option to estimate different linear parameters for each field
( \code{ collapseFixedEffect = FALSE } ) and in this case a table is not created because there is more than one estimate. See (\code{Omega} and \code{fixedEffectsCov}) in the \code{mKrig} object to build the standard errors.

\strong{plot} method  provides a panel of 4 diagnositic plots of the fit. 
Use \code{set.panel(2,2)} to see all 4 at once. The third plot gives the likelihood and
GCV functions as a function of lambda evaluated at the global MLE for theta.
This is based on the gird evaluations in the component MLEInfo$MLEProfileLambda.
The fourth
plot is a profile likelihood trace for theta  having maximized over lambda and is based on the component MLEInfo$MLEGrid. 

\strong{print} method prints the  \code{spatialProcessSummary} object of the fit, adding some details and explanations.

%%  ~~ If necessary, more details than the description above ~~
}
\value{
	An object of classes \code{mKrig} and  \code{SpatialProcess}. The main difference
	from mKrig is an extra component, \code{MLEInfo} that has the results of the profile likelihood grid
	evaluation over theta ( having maximizing lamdba), joint maximization over theta and lambda,
	and a grid evaluation over lambda with theta fixed at its MLE. 
}

\author{
Doug Nychka%%  ~~who you are~~
}

\seealso{
Tps, MLESpatialProcess, mKrigMLEGrid, mKrigMLEJoint, plot.Krig, predict.mKrig, predictSE.mKrig 
}
\examples{
data( ozone2)
# x is a two column matrix where each row is a location in lon/lat 
# coordinates
  x<- ozone2$lon.lat
# y is a vector of ozone measurements at day 16 a the locations. 
  y<- ozone2$y[16,]
  obj<- spatialProcess( x, y)
# summary of model
  summary( obj)
# diagnostic plots
  set.panel(2,2)
  plot(obj)
# plot 1 data vs. predicted values
# plot 2 residuals vs. predicted
# plot 3 criteria to select the smoothing
#        parameter lambda = sigma^2 / rho
#        the x axis has transformed lambda
#        in terms of effective degrees of freedom 
#        to make it easier to interpret
#        Note that here the GCV function is minimized
#        while the REML is maximzed. 
# plot 4 the log profile likelihood used to 
#        determine theta. 
#
# predictions on a grid 
  surface( obj)
#(see also predictSurface for more control on evaluation grid
# and plotting)
#  

\dontrun{
# working with covariates and filling in missing station data
# using an ensemble method
# see the example under  help(sim.spatialProcess) to see how to 
# handle a conditional simulation on a grid of predictions with 
# covariates. 
data(COmonthlyMet)
  fit1E<- spatialProcess(CO.loc,CO.tmin.MAM.climate, Z=CO.elev, 
                               theta.range= c(.25, 2.0) )
  set.panel( 2,2)                             
  plot( fit1E)
  
# conditional simulation at missing data
  notThere<- is.na(CO.tmin.MAM.climate )
  xp <- CO.loc[notThere,]
  Zp <- CO.elev[notThere]
  infill<- sim.spatialProcess( fit1E, xp=xp,
                      Z= Zp, M= 10)
#  
# interpretation is that these infilled values are all equally plausible 
# given the observations and also given the estimated covariance model
#  
# for extra credit one could now standardized the infilled values to have
# conditional mean and variance from the exact computations
#  e.g. predict( fit1E, xp=CO.loc[!good,],  Z= CO.elev[!good])
#  and  predictSE(fit1E, xp=CO.loc[!good,],  Z= CO.elev[!good])  
# with these standardization one would still preserve the correlations
# among the infilled values that is also important for considering them as a
# multivariate prediction.
# conditional simulation on a grid but not using the covariate of elevation
 fit2<- spatialProcess(CO.loc,CO.tmin.MAM.climate, 
                               theta.range= c(.25, 2.0) )
# note larger range parameter
# create 2500 grids using handy function
gridList <- fields.x.to.grid( fit2$x, nx=50,ny=50)
xGrid<- make.surface.grid( gridList)
ensemble<- sim.spatialProcess( fit2, xp=xGrid, M= 5)
# this is an "n^3" computation so increasing the grid size 
# can slow things down for computation 
image.plot( as.surface( xGrid, ensemble[1,]))
set.panel()
}

\dontrun{
data( ozone2)
# x is a two column matrix where each row is a location in lon/lat 
# coordinates
  x<- ozone2$lon.lat
# y is a vector of ozone measurements at day 16 a the locations. 
  y<- ozone2$y[16,]
# a comparison to using an exponential and Wendland covariance function
# and great circle distance -- just to make range easier to interpret.
    obj <- spatialProcess( x, y,
                              Distance = "rdist.earth")
	obj2<- spatialProcess( x, y,
	        cov.args = list(Covariance = "Exponential"), 
                              Distance = "rdist.earth" )
	obj3<- spatialProcess( x, y,
	        cov.args = list(Covariance = "Wendland",
	                        dimension  = 2,
	                                 k = 2),
	                          Distance = "rdist.earth")
# obj2 could be also be fit using the argument:
#   cov.args = list(Covariance = "Matern", smoothness=.5)
#	                          
# Note very different range parameters - BTW these are in miles
# but similar nugget variances. 
obj$pars
obj2$pars
obj3$pars
# since the exponential is Matern with smoothness == .5 the first two
# fits can be compared in terms of their likelihoods
# the REML value is slightly higher for obj verses obj2 (598.4  > 596.7)
# these are the _negative_ log  likelihoods so suggests a preference for the
# exponential model 
# 
# does it really matter in terms of spatial prediction?
set.panel( 3,1)
surface( obj)
US( add=TRUE)
title("Matern sm= 1.0")
surface( obj2)
US( add=TRUE)
title("Matern sm= .5")
surface( obj3)
US( add=TRUE)
title("Wendland k =2")
# prediction standard errors
# these take a while because prediction errors are based 
# directly on the Kriging weight matrix
# see mKrig for an alternative.
set.panel( 2,1)
out.p<- predictSurfaceSE( obj, nx=40,ny=40)
surface( out.p)
US( add=TRUE)
title("Matern sm= 1.0")
points( x, col="magenta")
#
out.p<- predictSurfaceSE( obj, nx=40,ny=40)
surface( out.p)
US( add=TRUE)
points( x, col="magenta")
title("Matern sm= .5")
}
set.panel(1,1)

\dontrun{
### comparison with GeoR
  data(ozone2)
  x<- ozone2$lon.lat
  y<- ozone2$y[16,]
  good<-!is.na(y)
  x1<- x[good,]
  y1<- y[good]
  
  obj<- spatialProcess( x, y, mKrig.args= list(m=1), smoothness = .5 )
  
  library( geoR)
  ml.n <- likfit(coords= x1, data=y1, ini = c(570, 3), nug = 50)
  # compare to 
  stuffFields<- obj$MLEInfo$MLEJoint$summary[c(1,3,4,5)]
  stuffGeoR<- c( ml.n$loglik, ml.n$phi, sqrt(ml.n$nugget),ml.n$sigmasq) 
  test.for.zero(  max(stuffFields/stuffGeoR), 1, tol=.004)    
}
}
\keyword{ spatial}
% __ONLY ONE__ keyword per line
