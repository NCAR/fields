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
\alias{spatialProcessSetDefaults}
\alias{plot.spatialProcess}
\alias{print.spatialProcess}
\alias{print.spatialProcessSummary}
\alias{summary.spatialProcess}
\alias{profileMLE}
\alias{confidenceIntervalMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Estimates a spatial process model. 
%%  ~~function to do ... ~~
}
\description{
	For a given covariance function  estimates the covariance parameters by maximum  maximum likelihood and then evaluates the 
	spatial model with these estimated parameters. The returned objet can be used for spatial prediction, conditional simulation, and profiling the likelihood function. 
	 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
spatialProcess(x, y,  weights = rep(1, nrow(x)),   Z = NULL,
                  mKrig.args = NULL,
                cov.function = NULL, 
                    cov.args = NULL,
                     parGrid = NULL, 
                      reltol = 1e-4,
                       na.rm = TRUE,
                     verbose = FALSE,
                        REML = FALSE, 
            cov.params.start = NULL,
                       gridN = 5,
               profileLambda = FALSE,
               profileARange = FALSE,
               profileGridN  = 15, 
                  gridARange = NULL,
                  gridLambda = NULL,
                 CILevel= .95,
                           ...)
\method{summary}{spatialProcess}(object,  ...)
\method{print}{spatialProcess}(x, digits = 4, ...)
\method{print}{spatialProcessSummary}(x, digits = 4, ...)
\method{plot}{spatialProcess}(x, digits = 4, which = 1:4, ...)

                           
spatialProcessSetDefaults( x, cov.function,
                                      cov.args,
                                      cov.params.start,
                                      mKrig.args,
                                      extraArgs=NULL,
                                      parGrid,
                                      gridN=5,
                                      verbose=FALSE)
    
confidenceIntervalMLE( obj, CILevel, verbose=FALSE)

profileMLE (obj, parName, parGrid=NULL, gridN=15,
                       cov.params.start=NULL, GCV=FALSE, REML=FALSE,
                       verbose=FALSE)


}

\arguments{
\item{x}{Observation locations}

\item{y}{Observation values}

\item{weights}{Weights for the error term (nugget) in units of
reciprocal variance.}

\item{Z}{A matrix of extra covariates for the fixed part of spatial
model.
E.g. elevation for fitting climate data over space. }

\item{CILevel}{Confidence level for intervals for the estimated
parameters.}

 \item{cov.args}{ A list specifying parameters and other components 
 of the
 covariance function. Default is not extra arguments required.]
 (But see the next item.). }
 
\item{cov.function}{A character string giving the name of the 
covariance function
  	 for the spatial component.
  	 If NULL, the default,  this is  filled in as 
  	 \code{stationary.cov} and then
  	 if \code{cov.args} is also NULL this is filled in as 
  	 \code{ list(Covariance = "Matern", smoothness = 1.0)} by the spatialProcessSetDefaults 
  	 function. 
  	 }
  	 
\item{cov.params.start}{ A list where the names are  parameter names
that appear in the covariance function. The values of each component
are assumed to be the starting values when optimizing to find MLEs. If
lambda does not appear
as additional
argument when calling spatialProcess it is added internally to this
list with the starting value .5.}

 \item{digits}{Number of significant digits in printed summary}

\item{extraArgs}{Extra arguments passed using the \dots device in R.
Typically these are extra covariance parameters specified in
\code{spatialProcess} and then passed to
\code{spatialProcessSetDefaults} }
 
\item{GCV}{A future argument not currently implemented. If TRUE will
find parameters by minimizing an approximate generalized 
cross-validation function. }
 
 \item{gridARange}{A grid for profiling  over the range parameter. If
 omitted, default is  based on a grid of profileGridN points centered at the
 MLE. 
}
\item{gridLambda}{A grid for profiling over lambda. }

 \item{gridN}{Number of grid points for intital fgrid search to find starting values. }
 
  	
\item{na.rm}{If TRUE NAs are removed from the data.}

\item{mKrig.args}{Arguments passed to the mKrig function.}

\item{obj}{A spatialProcess object returned from the spatialProcess function.}
\item{object}{ See \code{obj}.}

\item{parGrid}{ A data frame with the values of covariance parameters
to use as an initial grid search for starting values.}

\item{parName}{Text string that is the name of the parameter to 
 profile .}
 
 \item{profileARange}{If TRUE profile likelihood on aRange. Default is
TRUE if aRange is omitted.}

\item{profileGridN}{Number of grid points to use for profiling.}

\item{profileLambda}{If TRUE profile likelihood on lambda. This takes
extra time and is not necessary so the default is FALSE. }

 \item{reltol}{ Relative tolerance used in optim for convergence.}
 
 \item{REML}{ If TRUE the parameters are found by restricted maximum likelihood.}
 
 \item{verbose}{If TRUE print out intermediate information for debugging.}

  \item{\dots}{
  	Any other arguments that will be passed to the \code{mKrig} function and interpreted
  	as additional arguments to the covariance function.  This is a lazy way of specifying these. 
  	E.g. \code{aRange =.1} will set the covariance argument aRagne to .1. 
  }
  
  \item{which}{The vector 1:4 or any subset of 1:4, giving the plots
  to draw. See the description ofthese plots below.}

}
\details{
	This function makes many choices for the user in terms of defaults and it is 
	important to be aware of these. 
	The spatial model is 
	
	Y.k= P(x.k) + Z(x.k)\%*\%d2 + g(x.k) + e.k

where ".k" means subscripted by k, Y.k is the dependent variable
observed at location x.k. P is a low degree polynomial (default is a
linear function in the spatial coordinates) and Z is a matrix of covariates (optional) that
enter as a linear model the fixed part. g is a mean zero,
Gaussian stochastic process with a marginal variance of sigma and a
scale (or range) parameter, aRange. The measurement errors, e.k, are
assumed to be uncorrelated, normally distributed with mean zero and
standard deviation tau. If weights are supplied then the variance of e is assumed to be
\code{tau^2/ weights}.


Perhaps the most important aspect of this  function is that
 	the range (aRange), nugget (tau**2) and process variance (sigma) parameters
 	for the covariance are estimated by restricted maximum
 	likelihood and this is the model that is then  used for spatial
 	prediction. Geostatistics usaually refers to tau**2 + sigma as the 
 	"sill" and often these parameters are estimated by variogram fitting rather
 	than maximum likelihood. To be consistent with spline models and to focus
 	on the key part of model we reparametrize as lambda= tau**2/
 	sigma and sigma. Thinking about h as the spatial signal and e as the noise lambda
 	can be interpreted  as the noise to signal variance ratio in this spatial
 	context.(See the comparision with fitting the geoR model in the
examples section.)
 	
 	The likelihood and the cross valdiation function
 	can be concentrated to only depend on lambda and aRange and so
 	in reported the optimiztation of these two criterion we focus
 	on this form of the parameters.  Once lambda and aRange are
 	found, the MLE for sigma has a closed form and of course then
 	tau is then determined from lambda and sigma. 
 	
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
 	the user just wants to fix the range parameter aRange then
 	\code{Krig} can be used.
 	
 	NOTE: The defaults for the \code{optim} function used in MLESpatialProcess are:
 	
  \preformatted{
  list(method = "BFGS", 
       control=list(fnscale = -1,
                      ndeps = rep(log(1.1),length(cov.params.start)+1),
                     reltol = reltol,
                      maxit = 20))
}
                                     

There is always a hazard in providing a simple to use method that
makes many default choices for the spatial model. As in any analysis
be aware of these choices and try alternative models and parameter
values to assess the robustness of your conclusions. Also examine the
residuals to check the adequacy of the fit. See the examples below for
some help in how to do this easily in fields. Also see quilt.plot to
get an quick plot of a spatial field  to discern obvious spatial paterns.


\strong{summary} method forms a list of class \code{spatialProcessSummary} that has a
subset of information from the output object and also creates a table of the estimates
of the linear parameters in the fixed part of the model.
With replacated fields there is an option to estimate different linear parameters for each field
( \code{ collapseFixedEffect = FALSE } ) and in this case a table is not created because
there is more than one estimate. See (\code{Omega} and \code{fixedEffectsCov}) in the
\code{mKrig} object to build the standard errors.

\strong{plot} method  provides potentially  four diagnositic plots of the fit.Use the \code{which}
to pick and choose among them 
or use \code{set.panel} to see  them all.
The third and fourth plots, however, are only available if the profile computations been done. 
If lambda is profiled (\code{lambdaProfile} is not
\code{NULL} ) the third plot is the profile log
likelihood for lambda and with the GCV function on a
second verical scale.
This is based on the grid evaluations in the component 
\code{ lambdaProfile\$MLEProfileLambda} .
The fourth
plot is a profile log likelihood trace for aRange
based on \code{ aRangeProfile\$MLEProfileLambda}.


\strong{print} method prints the  \code{spatialProcessSummary} object of the fit, adding
some details and explanations.

\strong{spatialProcessSetDefaults} This is a useful way to fill in defaults for the function in one place. The main choices are choosing the Matern family, smoothness and a default fixed model (aka spatial drift). The grids for profiling are also created if they have not been supplied. 

%%  ~~ If necessary, more details than the description above ~~
}
\value{
	An object of classes \code{mKrig} and  \code{SpatialProcess}. The difference
	from mKrig are some extra  components. The more useful ones are listed belo
	
	
  \strong{MLESummary} A named array that has the fixed and estimated parameters along with likelihood values and some optim info.     
  
  \strong{profileSummaryLambda and profileSummaryARange} The output list from 
  mKrigMLEGrid for searching over over a grid of lambda and aRange. 
  
  \strong{CITable} Approximate confidence intervals based on the inverse hessian of the log likelihood function. 
  
  \strong{MLEInfo} A list that has a full documentation of the maximization including all parameters  and likelihood values that were tried by the optim function. 
  
   \strong{InitialGridSearch} Results from initial grid search to get good starting values for lambda and/or aRange.
        	
}

\author{
Doug Nychka%%  ~~who you are~~
}

\seealso{
\link{Tps},   \link{mKrigMLEGrid},
\link{mKrigMLEJoint},  \link{plot.Krig},  \link{predict.mKrig},
\link{predictSE.mKrig} 
}
\examples{
data( ozone2)
# x is a two column matrix where each row is a location in lon/lat 
# coordinates
  x<- ozone2$lon.lat
# y is a vector of ozone measurements at day 16. Note some missing values. 
  y<- ozone2$y[16,]
  
# artifically reduce size of data for a quick example to pass CRAN ...
  x<- x[1:75,]
  y<- y[1:75]
  
# lots of default choices  made here -- see gridN to increase 
# the number of points in grid searches for MLEs
# without specifying  lambda or aRange both are found in a robust 
# way uses grid searches 
# profiling over lambda and aRange  is not reuqired but completes the full 
# example.  Omit this for a faster computation.

  obj<- spatialProcess( x, y, profileLambda=TRUE, profileARange=TRUE)
  
# summary of model
  summary( obj)
  
# diagnostic plots
  set.panel(2,2)
  plot(obj)
  
# plot 1 data vs. predicted values
# plot 2 residuals vs. predicted
# plot 3 criteria to select the smoothing
#        parameter lambda = tau^2 / sigma
#        the x axis has log10 lambda
#        Note that here the GCV function is minimized
#        while the log profile likelihood  is maximzed. 
# plot 4 the log profile likelihood used to 
#        determine range parameter aRange. 
#

set.panel()
# predictions on a grid 
surface( obj, xlab="longitude", ylab="latitude")
US( add=TRUE, col="grey", lwd=2)
title("Predicted ozone (in PPB)  June 18, 1987 ")  
#(see also predictSurface for more control on evaluation grid, predicting
# outside convex hull of the data. and plotting)

# prediction standard errors, note two steps now to generate
# and then plot surface
look<- predictSurfaceSE( obj)
surface( look, xlab="longitude", ylab="latitude")
points( x, col="magenta")
title("prediction standard errors (PPB)")

# here  is a sanity check -- call spatialProcess with the MLEs found
# above, better get the same predictions!

objTest<- spatialProcess( x, y,
                          lambda=obj$MLESummary["lambda"],
                          aRange=obj$MLESummary["aRange"]
                          )
                          
test.for.zero(objTest$fitted.values, obj$fitted.values, 
            tag="sanity check" )

\dontrun{
# working with covariates and filling in missing station data
# using an ensemble method
# see the example under  help(sim.spatialProcess) to see how to 
# handle a conditional simulation on a grid of predictions with 
# covariates. 
data(COmonthlyMet)
  fit1E<- spatialProcess(CO.loc,CO.tmin.MAM.climate, Z=CO.elev, 
                               gridARange= seq(.25, 2.0, length.out=10)
			       )
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
# EXTRA CREDIT:  standardize the infilled values to have
# conditional mean and variance from the exact computations
#  e.g. predict( fit1E, xp=CO.loc[!good,],  Z= CO.elev[!good])
#  and  predictSE(fit1E, xp=CO.loc[!good,],  Z= CO.elev[!good])  
# with these standardization one would still preserve the correlations
# among the infilled values that is also important for considering them as a
# multivariate prediction.
# conditional simulation on a grid but not using the covariate of elevation
 fit2<- spatialProcess(CO.loc,CO.tmin.MAM.climate,
                        gridARange= seq(.25, 2.0, length.out=10)
                       )
# note larger range parameter
# create 2500 grid points using a handy fields function
gridList <- fields.x.to.grid( fit2$x, nx=50,ny=50)
xGrid<- make.surface.grid( gridList)
ensemble<- sim.spatialProcess( fit2, xp=xGrid, M= 5)
# this is an "n^3" computation so increasing the grid size 
# can slow things down for computation 
image.plot( as.surface( xGrid, ensemble[1,]))
set.panel()
}

\dontrun{
## changing the covariance model.

data(ozone2)

  x<- ozone2$lon.lat
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
 rbind( Whittle= obj$summary,
 Exp= obj2$summary,
Wendland= obj3$summary
)

# since the exponential is Matern with smoothness == .5 the first two
# fits can be compared in terms of their likelihoods
# the ln likelihood   value is slightly higher for obj verses obj2 (-613.9 >  -614.9)
# these are the _negative_ log  likelihoods so suggests a preference for the
# smoothness = 1.0 (Whittle)  model 
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
set.panel(1,1)
}


\dontrun{
### comparison with GeoR
data(ozone2)
    x<- ozone2$lon.lat
    y<- ozone2$y[16,]
    good<-!is.na(y)
    x1<- x[good,]
    y1<- y[good]
    
# spatialProcess takes longer because of grid search on aRange.    
    obj<- spatialProcess( x1, y1,
                          mKrig.args = list(m=1),
  			                  smoothness = .5)
    library( geoR)
    ml.n <- likfit(coords= x1, data=y1, ini = c(570, 3), nug = 50)
    # compare to 
    stuffFields<- obj$MLESummary[c("lnProfileLike.FULL",
                                   "aRange","tau","sigma2")]
    stuffGeoR<- c( ml.n$loglik, ml.n$phi, sqrt(ml.n$nugget),
                   ml.n$sigmasq) 
    test.for.zero(  max(stuffFields/stuffGeoR), 1, tol=.004)
   
}
}
\keyword{ spatial}
% __ONLY ONE__ keyword per line
