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
spatialProcessOLD <- function(x, y,  weights = rep(1, nrow(x)),   Z = NULL,
                  mKrig.args = NULL,
                cov.function = NULL, 
                  	cov.args = NULL,
                   gridARange = NULL,
                   gridLambda = NULL,
                      reltol = 1e-7,
                       na.rm = TRUE,
                  	 verbose = FALSE,
                        REML = FALSE, 
             confidenceLevel = .95,
            cov.params.start = NULL,
                       gridN = 5,
               profileLambda = FALSE,
               profileARange = FALSE,
               gridSearch= TRUE,
               parGrid= NULL, 
                           ...) {
 
  
# NOTE all ... information is assumed to be for the cov.args list
  GCV<- FALSE # top level placeholder to add GCV search capability once
# algorithm is stable 
# this default choice is used in the next level functions.
  
# set defaults based on the model passed. 
# obj is the output list and will be added to throughout the computation
   extraArgs<- list(...)
   
   if( REML&GCV){
     stop("Cannot optimize for both REML and GCV!")
   }
    
   obj<- spatialProcessSetDefaults(cov.function=cov.function, 
                                   cov.args=cov.args,
                                   cov.params.start=cov.params.start,
                                   mKrig.args = mKrig.args,
                                   extraArgs = extraArgs,
                                   profileLambda=profileLambda,
                                   profileARange=profileARange,
                                   gridARange=gridARange,
                                   gridLambda=gridLambda,
                                   gridN=gridN
                                   )
   #
   # obj$CASE 
   # 0 evaluate on passed cov parameters but MLEs for sigma, tau found from
   #   lambda
   #
   # 1 optimize loglikelihood over parameter lambda or  lambda and theta if
   #   not specified 
   #
   # 2 optimize loglikelihood over any parameters speficied in 
   #   cov.params.start but not in cov.args
   #
   # 3 grid search/profile over  a grid of aRange or use 
   #   gridARange as the grid
   #
   # 4 grid search/profile over  a grid of lambda
   #
   # 5 profile over aRange followed by profile over lambda.
   #
# Being in the cov.params.start list means a parameter will be optimized
# by maximum likelhood. If a  parameter is in the cov.args list it will be fixed in its value
# 
  if( verbose){
    cat(" The CASE:", obj$CASE, fill=TRUE)
    cat("Complete list of components in cov.args: ", "\n",
        names(obj$cov.args),fill=TRUE )
  }
####################################################################  
# CASE 0 small grid search over aRange and lambda
####################################################################  
   if( !is.null(parGrid)){
      cov.params.startTemp <- NULL 
     theGridSearch<- mKrigMLEGrid(x, y,  
                                  weights = weights, Z= Z, 
                                  mKrig.args = mKrig.args,
                                  cov.function = obj$cov.function, 
                                  cov.args  = obj$cov.args,
                                  par.grid = parGrid, 
                                  na.rm = na.rm,
                                  verbose = FALSE,
                                  REML = REML,
                                  GCV = GCV,
                                  reltol = reltol,
                                  cov.params.start = NULL)
     return(theGridSearch)
   }

####################################################################  
# CASE 3 and part I of CASE 5
####################################################################  
   
  if( (obj$CASE == 3) | (obj$CASE == 5) ){
    
   
    cov.params.startTemp <-  obj$cov.params.start 
    cov.params.startTemp$aRange<- NULL
    
    aRangeProfile<- mKrigMLEGrid(x, y,  
                                 weights = weights, Z= Z, 
                                 mKrig.args = mKrig.args,
                                 cov.function = obj$cov.function, 
                                 cov.args  = obj$cov.args,
                                 par.grid = list( aRange = obj$gridARange), 
                                 na.rm = na.rm,
                                 verbose = FALSE,
                                 REML = REML,
                                  GCV = GCV,
                                  reltol = reltol,
                     cov.params.start = cov.params.startTemp)
    indMax<- aRangeProfile$indMax
    aRangeMLEGrid<- aRangeProfile$summary[indMax,"aRange"]
    # refine starting values based on the grid search over aRange
    obj$cov.params.start$aRange<- aRangeMLEGrid
    lambdaMLE<- aRangeProfile$summary[indMax,"lambda"]
    aRangeMLE0<- aRangeProfile$summary[indMax,"aRange"]
    if( verbose){
      cat(lambdaMLE,aRangeMLE0, fill=TRUE )
      cat("ARangeProfile summary", fill=TRUE)
      print( aRangeProfile$summary )
    }
  }
   
####################################################################  
# profile/ grid search on lambda and PART II CASE 5
#################################################################### 
   
   if(  profileLambda | obj$CASE==4){
     lambdaCenter<- obj$cov.params.start$lambda
     if( is.null( gridLambda)){
       gridLambda = 10**( seq( -2,2, length.out=gridN))*lambdaCenter
     }
     cov.params.startTemp<- obj$cov.params.start
     cov.params.startTemp$lambda<- NULL
     if(verbose){
       cat("call to profile grid cov.args")
       print( obj$cov.args)
     }
     lambdaProfile<- mKrigMLEGrid(x, y,  
                                weights = weights, Z= Z, 
                             mKrig.args = mKrig.args,
                                cov.function = obj$cov.function, 
                              cov.args  = obj$cov.args,
                               par.grid = list( lambda = gridLambda), 
                                  na.rm = na.rm,
                                verbose = verbose,
                                   REML = REML,
                                    GCV = GCV,
                                    reltol= reltol,
                       cov.params.start = cov.params.startTemp
                       )
     indMax<- lambdaProfile$indMax
     lambdaMLEGrid<- lambdaProfile$summary[indMax,"lambda"]
     # refine starting values based on the grid search over aRange
     obj$cov.params.start$lambda<- lambdaMLEGrid
     if(verbose){
       cat("profile summary", fill=TRUE)
       print( lambdaProfile$summary)
     }
   } 
   
####################################################################  
# CASES 1 , 2 , 3, 4
#################################################################### 
   if(obj$CASE !=0 ){
     # optimze over all parameters 
     # where starting values are given or if
     # values in cov.args are omitted. 
     if(verbose){
       cat(" ********** starts used for final optimization", fill=TRUE)
       print( obj$cov.params.start)
     }
     MLEInfo <-mKrigMLEJoint(x, y,  weights = weights, Z= Z, 
                             mKrig.args = obj$mKrig.args,
                             cov.function = obj$cov.function, 
                             cov.args  = obj$cov.args,
                             na.rm = na.rm,
                             verbose = verbose,
                             reltol= reltol,
                             cov.params.start = obj$cov.params.start,
                             REML = REML,
                             GCV = GCV) 
   }
   
################################################################################
# final fit 
# now fit spatial model with MLE(s) 
# or the value(s) supplied in the call
# reestimate the other parameters for simplicity to get the
# complete mKrig object
####################################################################
   
# if all parameters are fixed -- don't mess with cov.args   
   if( obj$CASE == 0){
     obj$cov.argsFull <-  obj$cov.args
     
   }
   else{
      dupParameters<- match( names(MLEInfo$pars.MLE ), names(cov.args) )
     if( all( is.na(dupParameters)) ){
       obj$cov.argsFull<- c( obj$cov.args,
                             as.list(MLEInfo$pars.MLE) )
     }
   }
  if( verbose){
    cat( "Names cov.args:","\n", names( obj$cov.args), fill=TRUE)
    cat( "Names cov.params.start:","\n", names( obj$cov.params.start),fill=TRUE)
    cat( "Names argsFull:","\n", names( obj$cov.argsFull),fill=TRUE  )
  }
  mKrigObj <- do.call( "mKrig", 
	                c( list(x=x,
	                        y=y,
	                  weights=weights,
	                        Z=Z),
	                  obj$mKrig.args,
	             list( na.rm=na.rm),
	             list(cov.function = obj$cov.function),
	                  obj$cov.argsFull 
	            	)
             	)
####################################################################
# sort out output object based on the different cases
####################################################################
  if( obj$CASE==0){
    MLEInfo <- NULL
    aRange.CI<- NA
    lambda.CI<- NA
    aRangeProfile<- NULL
    lambdaProfile<- NULL
    MLESummary <- mKrigObj$summary
  }
  
  if( obj$CASE==1 | obj$CASE==2){
    aRangeProfile<- NULL
    lambdaProfile<- NULL
    aRange.CI<- NA
    lambda.CI<- NA
    MLESummary <- MLEInfo$summary
  }
  if( obj$CASE==3){
   aRange.CI<- profileCI( aRangeProfile, "aRange", confidenceLevel)
   lambdaProfile<- NULL
   lambda.CI<- NA
   MLESummary <- MLEInfo$summary
  }
  if( obj$CASE==4){
    aRangeProfile<- NULL
    aRange.CI<- NA
    lambda.CI<- profileCI( lambdaProfile, "lambda", confidenceLevel)
    MLESummary <- MLEInfo$summary
  }
  if( obj$CASE==5){
    aRange.CI<- profileCI( aRangeProfile, "aRange", confidenceLevel)
    lambda.CI<- profileCI( lambdaProfile, "lambda", confidenceLevel)
    MLESummary <- MLEInfo$summary
  }
# combine everything into the output list  
	obj <- c(obj, mKrigObj, 
	                 list(MLESummary=MLESummary,
	                      MLEInfo = MLEInfo,
	                      aRangeProfile = aRangeProfile, 
	                      lambdaProfile = lambdaProfile,
	                      lambda.CI =  lambda.CI,
	                      aRange.CI =  aRange.CI
	                      )
	        )
# replace call in mKrig  object with the top level one
# from spatialProcess
  obj$call<- match.call()	
	class(obj) <- c( "spatialProcess","mKrig")
 
	return(obj)
}
