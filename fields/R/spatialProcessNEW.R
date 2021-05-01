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
spatialProcessNEW <- function(x, y,  weights = rep(1, nrow(x)),   Z = NULL,
                  mKrig.args = NULL,
                cov.function = NULL, 
                  	cov.args = NULL,
                      parGrid = NULL, 
                      abstol = 1e-4,
                       na.rm = TRUE,
                  	 verbose = FALSE,
                        REML = FALSE, 
             confidenceLevel = .95,
            cov.params.start = NULL,
                       gridN = 5,
               profileLambda = FALSE,
               profileARange = FALSE,
                  gridARange = NULL,
                  gridLambda = NULL,
                doGridSearch = TRUE,
            
                           ...) {
 
#  THE RULES: 
# Being in the cov.params.start list means a parameter will be optimized
# by maximum likelhood. 
# If a  parameter is in the cov.args list it will be fixed in its value
#  
  
# NOTE all the  ... extra arguents are assumed to be for the cov.args list
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
                                   parGrid=parGrid,
                                   gridN = gridN,
                                   doGridSearch=doGridSearch,
                                   verbose=verbose)
   #
   # obj$CASE 
   # 0 evaluate on passed cov parameters but MLEs for sigma, tau found from
   #   lambda
   #
   # 1 optimize loglikelihood over any parameters specified in 
   #   cov.params.start but not in cov.args
   #
   # 2  grid search over parameters using parGrid and generating starting values for 
   #   for the MLEs 
   #
   # 3 a default grid search over lambda and possibly aRange if
   #   these are not specified. 
   #    This is just a convenient special case 
   # 4 profile over  lambda and/or aRange 
   # this is fairly computationally demanding. 
   
  if( verbose){
    cat(" The CASE:", obj$CASE, fill=TRUE)
    cat("Complete list of components in cov.args: ", "\n",
        names(obj$cov.args),fill=TRUE )
  }
   
   
####################################################################  
# CASE 2 grid search for starting values 
####################################################################  
   
  if( (obj$CASE == 2 )  ){
    InitialGridSearch<- mKrigMLEGrid(x, y,  
                             weights = weights,
                                   Z = Z, 
                          mKrig.args = mKrig.args,
                        cov.function = obj$cov.function, 
                           cov.args  = obj$cov.args,
                            par.grid = parGrid, 
                               na.rm = na.rm,
                             # verbose = verbose,
                                REML = REML,
                                 GCV = GCV,
                    cov.params.start = cov.params.start)
  # use grid search to set starting values
    parN<- names( parGrid)
    print( cov.params.start)
    cov.params.start[parN] <- parGrid[InitialGridSearch$indMax, parN]
    print( cov.params.start)
  }
   else{
     InitialGridSearch<-NULL
   }
   
####################################################################  
# profile on lambda and/or theta
#################################################################### 
   
   if(  profileLambda){
    print("placeholder profileLambda")
   } 
   
   if(  profileARange){
     print("placeholder profileLambda")
   } 
   
####################################################################  
# CASES 1 , 2 , 3, 4
#################################################################### 
   if(obj$CASE !=0 ){
     # optimze over all parameters 
     # where starting values are given or if
     # values in cov.args are omitted. 
     MLEInfo <-mKrigMLEJoint(x, y,  weights = weights, Z= Z, 
                             mKrig.args = obj$mKrig.args,
                             cov.function = obj$cov.function, 
                             cov.args  = obj$cov.args,
                             na.rm = na.rm,
                             verbose = verbose,
                             cov.params.start = cov.params.start,
                             REML = REML,
                             GCV = GCV,
                             hessian = TRUE) 
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
    cat( "Names cov.params.start:","\n", names(cov.params.start), fill=TRUE)
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
  
# combine everything into the output list  
	obj <- c(obj, mKrigObj, 
	                 list(MLESummary=MLESummary,
	                      MLEInfo = MLEInfo,
	                      InitialGridSearch = InitialGridSearch,
	                      aRangeProfile = aRangeProfile, 
	                      lambdaProfile = lambdaProfile
	                      )
	        )
# replace call in mKrig  object with the top level one
# from spatialProcess
  obj$call<- match.call()	
	class(obj) <- c( "spatialProcess","mKrig")
 
	return(obj)
}
