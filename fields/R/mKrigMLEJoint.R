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
                          cov.function = "stationary.cov",
                              cov.args = NULL, 
                      cov.params.start = NULL,
                            optim.args = NULL,
                                reltol = 1e-6,
                          parTransform = NULL,
                                  REML = FALSE, 
                                  GCV  = FALSE,
                               hessian = FALSE, 
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
 
# precompute distance matrix if possible so it only needs to be computed once
  supportsDistMat = supportsArg(cov.function, "distMat")
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
# if lambda is  then it has been added to mKrig.args 
# if lambda.start then it is part of the parameter names and will 
# added in the cov.args list 
  mKrig.args <- c(list(x = x, y = y, weights = weights, Z = Z),
                   mKrig.args,
                  list(cov.function=cov.function) 
                  )
  if( verbose){
    cat("Info from mKrigJoint:",fill=TRUE)
    cat("Argument names in full mKrig.args: ",  fill=TRUE)
    print( names(mKrig.args) )
    cat("Full cov.args names:\n ", names( cov.args), fill=TRUE)
    cat("Parameters to optimze: ", parNames, length( parNames), fill=TRUE)
    cat("Starting values  (cov.params.start) :  ",  fill=TRUE)
    print( cov.params.start)
  }
  
  if(!GCV){
    nameCriterion<- ifelse( !REML,
                          "lnProfileLike.FULL",
                          "lnProfileREML.FULL" )
  }
  else{
    nameCriterion<-"GCV"
  }
  if(verbose){
    cat("nameCriterion", nameCriterion, fill=TRUE)
  }
  callOptim<-  !is.null(cov.params.start) & length(parNames) > 0
  
##########################################################################
###### main if else block
##########################################################################  
  if(callOptim){
   
    #########################################################
    ###  actually do something 
    #set default optim.args if necessary
    # abstol is anticipating this is a log likelihood so differencs of 1e-4 are not appreciable
    # 
    # setup control parameters for optim ....
    if( length(cov.params.start)==0){
      stop("On no! Found zero parameters to optimize!")
    }
    if(is.null(optim.args)){
      # number of step sizes either include lambda as a parameter or not.   
      ndeps<- rep(log(1.1), length( parNames) ) 
      optim.args = list(method = "BFGS", hessian = hessian,
                        control=list(fnscale = -1,
                                     ndeps = ndeps, 
                                     reltol = reltol,
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
  capture.evaluations <- NULL
 
  capture.env <- environment()
# call to optim with initial start (default is log scaling )
# 
  init.start <- parTransform( unlist(c(cov.params.start)), inv=FALSE)
  
  if( verbose){
  cat("Transformed starting values ","\n", init.start, fill=TRUE)
  }
# and now the heavy lifting ...
# optimize over (some) covariance parameters and possibly lambda
  optimResults <- do.call(optim, c(list(par=init.start),
                            list(mKrigJointTemp.fn),
                                         optim.args,
                           list(  parNames = parNames,
                              parTransform = parTransform,
                                mKrig.args = mKrig.args,
                                  cov.args = cov.args, 
                               capture.env = capture.env,
                                      REML = REML,
                                      GCV  = GCV)
                           )
                  )
  # reformat the  optim results
  lnLikeEvaluations <- capture.evaluations[-1,] # first row is just NAs
  if( verbose){
    cat("Captured  evaluations from optim: ", fill=TRUE)
    print(lnLikeEvaluations)
  }
  #ind<- which(lnLikeEvaluations[ , nameCriterion]
  #                         ==  optimResults$value )
  #ind<- max(ind)
  optim.counts <- optimResults$counts
  parOptimum<-    parTransform(optimResults$par, inv=TRUE)
  names( parOptimum)<- parNames
# update parameter start values with converged ones
# note that lambda may be part of this 
  cov.params.final<- as.list(parOptimum)
}
##########################################################################
###### end if block
########################################################################## 
else{
  # no optimization required, just setup for a single evaluation
  cov.params.final <- NULL
  lnLike.eval      <-  NA
  optimResults     <- NULL
  optim.counts     <- NA
  parOptimum       <- NULL
  lnLikeEvaluations<- NULL
  if( verbose){
    cat("fast return in mKrigMLEJoint",  fill=TRUE)
  }
}
#########################################################
### just evaluate
### at final parameters and also find the trace and GCV 
#########################################################
  cov.args.final<- c( cov.args, cov.params.final)
# mKrig.args$find.trA <- TRUE
  fastObject   <- do.call("mKrig",
                          c(mKrig.args, 
                          cov.args.final) )$summary
  
######################################################### 
  summary <- c( fastObject, optim.counts,
                parOptimum
                )

  out =list( summary = summary, 
            pars.MLE = parOptimum ,
        parTransform = parTransform,
        optimResults = optimResults , 
   lnLikeEvaluations = lnLikeEvaluations)
           
  return(out)
}



# Define the objective function as a tricksy call to mKrig
# if y is a matrix of replicated data sets use the log likelihood for the complete data sets
 mKrigJointTemp.fn <- function(parameters,
                              mKrig.args, cov.args, parTransform, parNames,
                              REML=FALSE,
                              GCV= FALSE,
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
  #
  # Don't find approx trace of A for log likelihood evaluations. 
  mKrig.argsTemp<- mKrig.args
  if( !GCV){
    mKrig.argsTemp$find.trA<- FALSE 
  }
  #  the summary vector from the mKrig object has everything we need ...
  hold <- do.call("mKrig", c(mKrig.argsTemp,cov.args.temp) )$summary
  # add this evalution to an object (i.e. here a matrix) in the calling frame
  temp.eval <- get("capture.evaluations", envir=capture.env)
  assign("capture.evaluations", 
         rbind(temp.eval,
            c(parTransform(parameters, inv=TRUE), hold)
            ), 
         envir = capture.env)
  
  if( !GCV){
    objectiveFunction <- ifelse(REML, 
                               hold["lnProfileREML.FULL"],
                               hold["lnProfileLike.FULL"])
   }
    else{
      objectiveFunction <- -1* hold["GCV"]
    }
 # cat("objective Function",objectiveFunction, fill=TRUE)
   
    return( objectiveFunction)
}

