# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2018
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
spatialProcess <- function(x, y,  weights = rep(1, nrow(x)),   Z = NULL,
                  mKrig.args = list( m=2),
                cov.function = NULL, 
                  	cov.args = NULL,
                   gridTheta = NULL, 
                      abstol = 1e-4,
                       na.rm = TRUE,
                  	 verbose = FALSE,
                        REML = FALSE, 
             confidenceLevel = .95,
            cov.params.start = NULL,
                       gridN = 10,
                           ...) {
# NOTE all ... information is assumed to be for the cov.args list
# overwrite the default choices (some R arcania!)
# default covariance is  Matern  with smoothness 1.
# linear regression model also added as the default. 
  
  ind<- match( names( cov.args), names(list(...) ) )
  cov.args <- c(cov.args[is.na(ind)], list(...))
  
  if( is.null( cov.function)){
    cov.function <- 'stationary.cov'
    if( is.null(cov.args) ){
      cov.args<- list()
    }
    if( is.null(cov.args$Covariance )){
        cov.args$Covariance<- "Matern"
        if( is.null(cov.args$smoothness )){
          cov.args$smoothness<- 1.0
        }
    }
  }  
#
# theta and lambda  are  handled specially because are almost always 
# estimated and this will simplify the call in  this top level function 
# NOTE that being in the cov.params.start list means a parameter will be optimized
# by maximum likelhood
#
# If a  parameter is in the cov.args list it will be fixed in its value
# 
  if( is.null( cov.args$lambda) & is.null(cov.params.start$lambda) ){
    cov.params.start$lambda <- .5
  }
  
  if( is.null( cov.args$theta) & is.null(cov.params.start$theta) ){
    cov.params.start$theta <- NA
  }
  
  if( !is.null( cov.params.start$theta ) & is.null( cov.params.start$lambda ) ){
    stop(" optimizing over theta for fixed lambda not supported")
  }
  
  if( verbose){
    cat("extra arguments from ... " , names( list(...)) , fill=TRUE)
    cat(" full list from cov.args: ", names(cov.args) )
  }
# NOTE: if theta is specified then just do a simple optimization
# if not  ... first a grid search in log theta followed by an full optimization 
  if( !is.null(cov.args$theta)){
###########################################################
#   case to   optimze over all parameters with theta fixed 
# parameter 
# list to optimize over if this list is NULL then function 
# will just evaluate at the parameters 
    cov.params.startTemp <-  cov.params.start
    MLEInfo <-mKrigMLEJoint(x, y,  weights = weights, Z= Z, 
              mKrig.args = mKrig.args,
                 cov.fun = cov.function, 
               cov.args  = cov.args,
                   na.rm = na.rm,
                 verbose = verbose,
        cov.params.start = cov.params.startTemp,
                    REML = REML) 
    if( is.null( cov.params.startTemp$lambda)){
      lambdaModel <- cov.args$lambda
      lambda.MLE <- NA
     }
    else{
    lambdaModel <- MLEInfo$pars.MLE["lambda"] 
     lambda.MLE <- lambdaModel
      theta.MLE <- NA
    theta.CI    <- NA
    thetaModel  <- cov.args$theta
  MLEGridSearch <- NULL
    }
  }
  else{
#  
	MLEGridSearch <- MLESpatialProcess(x, y, weights = weights, Z=Z, 
	                                mKrig.args = mKrig.args,
	                              cov.function = cov.function, 
	                                  cov.args = cov.args,
	                                 gridTheta = gridTheta,
	                                   	 gridN = gridN,
	                                    abstol = abstol,
                                  		verbose = verbose,
	                                       REML = REML,
	                           cov.params.start = cov.params.start
	                                       	 )
	
	lambda.MLE <- MLEGridSearch$MLEJoint$pars.MLE["lambda"] 
	lambdaModel<- lambda.MLE
	theta.MLE  <-   MLEGridSearch$MLEJoint$pars.MLE["theta"]
	thetaModel <-  theta.MLE
# approximate confidence interval for theta 
	thetaGrid<-     MLEGridSearch$MLEGrid$par.grid$theta
  lgProfileLike<- MLEGridSearch$MLEGrid$summary[,"lnProfileLike.FULL"]
	splineFit<- splint(thetaGrid, lgProfileLike, nx=500)
	cutLevel<- max(splineFit$y ) - qchisq(confidenceLevel, 1) / 2
	ind<- splineFit$y> cutLevel
	lower<- min( splineFit$x[ ind] )
	upper<- max(splineFit$x[ ind])
	theta.CI = c( lower, upper)
	MLEInfo<- MLEGridSearch$MLEJoint
  }
#  
	if( verbose){
	  cat("Summary from joint optimization", fill=TRUE)
	  print( MLEInfo$summary )
	  print( MLEInfo$pars.MLE)
	  print(theta.CI )
	}
  
################################################################################
# final fit 
# now fit spatial model with MLE for theta (range parameter)
# or the value supplied in the call and any other parameters that have been 
# been estimated. 
# re-estimate the other parameters for simplicity to get the
# complete mKrig object
  cov.argsFull <-  cov.args
  dupParameters<- match( names(MLEInfo$pars.MLE ), names(cov.args) )
  if( all( is.na(dupParameters)) ){
  cov.argsFull<- c( cov.argsFull, as.list(MLEInfo$pars.MLE) )
  }
  
  if( verbose){
    print( names( cov.args))
    print( names( cov.params.start))
    print( names(cov.argsFull) )
  }
  
  obj <- do.call( "mKrig", 
	                c( list(x=x,
	                        y=y,
	                  weights=weights,
	                        Z=Z),
	                  mKrig.args,
	          list( na.rm=na.rm),
	             list(cov.function = cov.function),
	                  cov.argsFull 
	            	)
	)
	obj <- c(obj,   
	                 list(   
	                       mKrig.args = mKrig.args,
	                     cov.function = cov.function, 
	                         cov.args = cov.args,
	                 cov.params.start = cov.params.start ),
	                 list(  
	                      MLEInfo = MLEInfo,
	                MLEGridSearch = MLEGridSearch,
	                   thetaModel = thetaModel,
	                    theta.MLE = theta.MLE,
	                     theta.CI = theta.CI,
	              confidenceLevel = confidenceLevel,
	                  lambdaModel = lambdaModel,
	                   lambda.MLE = lambda.MLE,
	                      summary = MLEInfo$summary)
	        )
# replace call to mKrig with this top level one
  obj$call<- match.call()	
	class(obj) <- c( "spatialProcess","mKrig")
 
	return(obj)
}
