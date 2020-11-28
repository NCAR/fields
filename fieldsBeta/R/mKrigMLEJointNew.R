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

mKrigMLEJoint<- function(x, y, weights = rep(1, nrow(x)),  Z = NULL,
                            mKrig.args = NULL,
                                 na.rm = TRUE,
                               cov.fun = "stationary.cov",
                              cov.args = NULL, 
                      cov.params.start = NULL,
                            optim.args = NULL,
                                abstol = 1e-4,
                          parTransform = NULL,
                              find.trA = TRUE, 
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
  # main way to keep track of parameters to optimize  
  # lambda is  included  if lambda.fixed is NULL
  # (most of the time one would want to optimize over lambda parameter)
    parNames <-   names(cov.params.start)
  if( verbose){
    cat("passed mKrig arguments: ", fill=TRUE)
    print( mKrig.args) 
    cat("parameters to optimze: ", parNames, fill=TRUE)
    cat("cov.params.start:  ",  fill=TRUE)
    print( cov.params.start)
  }
    #cat("starts", fill=TRUE)
    #print( cov.params.start)
# precompute distance matrix if possible so it only needs to be computed once
  supportsDistMat = supportsArg(cov.fun, "distMat")
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
  #add precomputed distance matrix to the  cov.args 
  cov.args = c(cov.args, list(distMat=distMat, onlyUpper=TRUE))
  }
# these are all the arguments needed to call mKrig except  cov.args
# if lambda.fixed then it has been added to mKrig.args 
# if lambda.start then it is part of the parameter names and will 
# added in the cov.args list 
  mKrig.args <- c(list(x = x, y = y, weights = weights, Z = Z),
                   mKrig.args,
                  list(cov.fun=cov.fun) 
                  )
  nameCriterion<- ifelse( !REML,
                          "lnProfileLike.FULL",
                          "lnProfileREML.FULL" )
##########################################################################
###### main if else block
##########################################################################  
  if(  !is.null(cov.params.start) ){
   
    #########################################################
    ###  actually do something 
    #set default optim.args if necessary
    # abstol is anticipating this is a log likelihood so differencs of 1e-4 are not appreciable
    # 
    # setup control parameters for optim ....
    if( length(parNames)==0){
      stop("Found zero parameters to optimize")
    }
    if(is.null(optim.args)){
      # number of step sizes either include lambda as a parameter or not.   
      ndeps<- rep(log(1.1), length( parNames) ) 
      optim.args = list(method = "BFGS", 
                        control=list(fnscale = -1,
                                     ndeps = ndeps, 
                                     abstol = abstol,
                                     maxit = 20)
      )
    }
    if( is.null(parTransform)){
      # parTransform: log/exp
      parTransform<- function( ptemp, inv=FALSE){
        if( !inv){ log( ptemp)}
        else{
          exp(ptemp)
        }
      }
    }     
  #
  capture.evaluations <- matrix(NA, ncol = length(parNames) + 4 , nrow = 1,
       dimnames = list(NULL,
                       c(   parNames,
                           "sigma.MLE",
                         "tau.MLE", 
                "lnProfileLike.FULL",
                "lnProfileREML.FULL")
                         )
                          ) 
  capture.env <- environment()
# call to optim with initial start (default is log scaling )
# 
  init.start <- parTransform( unlist(c(              cov.params.start)), inv=FALSE)
  if( verbose){
  cat("cov.params.start", init.start, fill=TRUE)
  print(cov.params.start )
  cat("initial parameter starts", init.start, fill=TRUE)
  }
# and now the heavy lifting ...
# optimize over (some) covariance parameters and possibly lambda
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
  # reformat the  optim results
  lnLike.eval <- capture.evaluations[-1,] # first row is just NAs
  ind<- which( lnLike.eval[ , nameCriterion]
                           == optimResults$value )
  ind<- max( ind)
  optim.counts <- optimResults$counts
  parOptimum<-    parTransform(optimResults$par, inv=TRUE)
  names( parOptimum)<- parNames
  optim.valueTest<-   optimResults$value
# update parameter start values with converged ones
# note that lambda may be part of this 
  cov.params.final<- as.list(parOptimum)
}
##########################################################################
###### end if block
########################################################################## 
else{
  # no optimization required, just setup for a single evaluation
  cov.params.final<- NULL
  lnLike.eval     <-  NA
  optimResults    <- NA
  parOptimum<- NULL
  optim.counts<- c( NA, NA)
  if( verbose){
    cat("fast return in mKrigMLEJoint",  fill=TRUE)
  }
}
#########################################################
### just evaluate
### at final parameters and also find the trace and GCV 
#########################################################
  cov.args.final<- c( cov.args, cov.params.final)
  fastObject   <- do.call("mKrig",
                        c(mKrig.args, list( find.trA= find.trA), cov.args.final) )
  tau.MLE    <- fastObject$tau.MLE.FULL
  sigma.MLE      <- fastObject$sigma.MLE.FULL
  optim.value  <- ifelse(!REML,  fastObject$lnProfileLike.FULL,
                         fastObject$lnProfileREML.FULL)
  eff.df <- fastObject$eff.df 
  GCV<- fastObject$GCV
######################################################### 
  summary <- c( optim.value, 
                parOptimum,
                tau.MLE,
                sigma.MLE,
                optim.counts,
                eff.df,
                GCV)
  names(summary) <- c(nameCriterion, parNames, 
                      "tauMLE", "sigmaMLE", "funEval", "gradEval","eff.df", "GCV")
  
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
  # NOTE  cov.args.temp can also include lambda as a component. 
  # due to the flexibility in arguments to mKrig  ( I.e. the ... argument)
  # any covariance arguments in cov.args.temp are matched to existing mKrig arguments
  # in particular lambda, the remaining unmatched arguments are assumed to be for the
  # covariance function and used within mKrig in the call to cov.fun
  
  hold <- do.call("mKrig", c(mKrig.args,
                             cov.args.temp))
  hold <-  hold[c("sigma.MLE.FULL",
                "tau.MLE.FULL",
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

