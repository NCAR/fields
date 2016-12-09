# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
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

mKrigMLEJoint <- function(x, y, weights = rep(1, nrow(x)),  Z = NULL,
                            mKrig.args = NULL,
                                 na.rm = TRUE,
                               cov.fun = "stationary.cov", cov.args=NULL, 
                          lambda.start = .5,
                      cov.params.start = NULL,
                            optim.args = NULL,
                                abstol = 1e-4,
                          parTransform = NULL, 
                                  REML = FALSE, 
                               verbose = FALSE) {
  # overwrite basic data to remove NAs this has be done in case distance 
  # matrices are precomputed (see below)
  if( na.rm){
    obj<- mKrigCheckXY(x, y, weights, Z, na.rm)
    x<- obj$x
    y<- obj$y
    weights<- obj$weights
    Z<- obj$Z
  }
  #set default optim.args if necessary
  # abstol is anticipating this is a likelihood so differencs of 1e-4 are not appreciable
  # 
  if(is.null(optim.args)){
    optim.args = list(method = "BFGS", 
                      control=list(fnscale = -1,
                                     ndeps = rep(log(1.1),length(cov.params.start)+1), 
                                    abstol = abstol,
                                     maxit = 20)
                      )
 }
# main way to keep track of parameters to optimize -- lambda always included  
parNames<- c( "lambda", names(cov.params.start))
if( is.null(parTransform)){
  # parTransform: log/exp
  parTransform<- function( ptemp, inv=FALSE){
    if( !inv){ log( ptemp)}
    else{
      exp(ptemp)
    }
  }
}
########bug
if(verbose){
  cat("parameters to optimze: ", parNames, fill=TRUE)
}
#check which optimization options the covariance function supports
  supportsDistMat = supportsArg(cov.fun, "distMat")
#precompute distance matrix if possible so it only needs to be computed once
  if(supportsDistMat & is.null( cov.args$distMat)) {
    #Get distance function and arguments if available
    #
    Dist.fun= c(cov.args)$Distance
    Dist.args=c(cov.args)$Dist.args
    
    #If user left all distance settings NULL, use rdist with compact option.
    #Use rdist function by default in general.
    #
    if(is.null(Dist.fun)) {
      Dist.fun = "rdist"
      if(is.null(Dist.args))
        Dist.args = list(compact=TRUE)
    }
  distMat = do.call(Dist.fun, c(list(x), Dist.args))
  #set cov.args for optimal performance
  cov.args = c(cov.args, list(distMat=distMat, onlyUpper=TRUE))
  }
# these are all the arguments needed to call mKrig except lambda and cov.args
  mKrig.args <- c(list(x = x, y = y, weights = weights, Z = Z),
                   mKrig.args,
                  list(cov.fun=cov.fun) 
                  )
# reset switch so trace is not found for each evaluation of the likelihood.   
  mKrig.args$find.trA = FALSE
# output matrix to summarize results
  ncolSummary = 7 + length(parNames)
  summary <- matrix(NA, nrow = 1, ncol = ncolSummary)
  dimnames(summary) <- list(NULL, c("EffDf", "lnProfLike", "GCV", "sigma.MLE", 
                                    "rho.MLE", parNames, 
                                    "counts eval","counts grad"))

  lnProfileLike.max <- -Inf
  
  #
  # optimize over (some) covariance parameters and lambda
  capture.evaluations <- matrix(NA, ncol = length(parNames) + 4 , nrow = 1,
       dimnames = list(NULL,
                       c(   parNames,
                           "rho.MLE",
                         "sigma.MLE", 
                "lnProfileLike.FULL",
                "lnProfileREML.FULL")
                         )
                          ) 
  capture.env <- environment()
# call to optim with initial start (default is log scaling )
  init.start <- parTransform( unlist(c(lambda.start, cov.params.start)), inv=FALSE)
#  cat("init.start",  init.start, fill=TRUE)
  optimResults <- do.call(optim, c(    list(par=init.start),
                            list(mKrigJointTemp.fn),
                                         optim.args,
                           list(  parNames = parNames,
                              parTransform = parTransform,
                                mKrig.args = mKrig.args,
                                  cov.args = cov.args, 
                               capture.env = capture.env,
                                      REML = REML)
                           )
                  )
#get optim results
  optim.counts <- optimResults$counts
  parOptimum<- parTransform(optimResults$par, inv=TRUE)
# first row is just NAs  
  lnLike.eval <- capture.evaluations[-1,]
 
  nameCriterion<- ifelse( !REML,
                          "lnProfileLike.FULL",
                          "lnProfileREML.FULL" )
  ind<- which( lnLike.eval[ , nameCriterion]
                       == optimResults$value )
  ind<- max( ind)
  # below is an aspect from optim I dont understand and thought to flag
  #if( length(ind)!=1 ){
  #     cat( "Weirdness in optimization. See lnLike.eval rows: ", ind,
  #         fill=TRUE )
  # ind<- max( ind)
  #}
# save results of the best covariance model evaluation in a neat table
  summary <- c(          optimResults$value, 
                                 parOptimum,
               lnLike.eval[ind,"sigma.MLE"],
                 lnLike.eval[ind,"rho.MLE"],
               optim.counts)
  names(summary) <- c(nameCriterion, parNames, 
                      "sigmaMLE", "rhoMLE", "funEval", "gradEval")
  out = c( list(summary=summary, lnLike.eval = lnLike.eval, optimResults=optimResults,
                    pars.MLE=parOptimum, parTransform = parTransform))
  return(out)
}

# Define the objective function as a tricksy call to mKrig
# if y is a matrix of replicated data sets use the log likelihood for the complete data sets
 mKrigJointTemp.fn <- function(parameters,
                              mKrig.args, cov.args, parTransform, parNames,
                              REML=FALSE,
                              capture.env) {
  # optimization is over a transformed scale ( so need to back transform for mKrig)
  tPars<- parTransform( parameters, inv=TRUE)
  names( tPars)<- parNames
  #get all this eval's covariance arguments using the input parameters
  cov.args.temp = c(cov.args, tPars)
  # NOTE: FULL refers to estimates collapsed across the replicates if Y is a matrix
  # assign to hold the last mKrig object
  hold <- do.call("mKrig", c(mKrig.args,
                             cov.args.temp))
  
  hold = hold[c("rho.MLE.FULL",
                "sigma.MLE.FULL",
                "lnProfileLike.FULL",
                "lnProfileREML.FULL"
                )]
  
  # add this evalution to an object (i.e. here a matrix) in the calling frame
  temp.eval <- get("capture.evaluations", envir=capture.env)
  assign("capture.evaluations", rbind(temp.eval,
                                      c(parTransform(parameters, inv=TRUE),
                                        unlist(hold))), 
         envir = capture.env)
  if( !REML){
  return(hold$lnProfileLike.FULL)
  }
  else{
  return(hold$lnProfileREML.FULL)
  }  
}

