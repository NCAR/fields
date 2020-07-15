
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

MLESpatialProcess <- function(x, y, weights = rep(1, nrow(x)), Z = NULL,
                            mKrig.args = NULL,
                            cov.function = "stationary.cov", 
                            cov.args = list(Covariance = "Matern",
                                            smoothness = 1), 
                             gridTheta = NULL,
                                 gridN = 20,
                            optim.args = NULL,
                                 na.rm = TRUE,
                               verbose = FALSE,
                               abstol  = 1e-4,
                                  REML = FALSE,
                      cov.params.start = NULL,
                               ...) {
  if( verbose){
    cat(" MLESpatialProcess extra arguments:" , fill=TRUE)
     print( names( list(...)))
  }
  # combine  list(...) with cov.args and omit duplicates but favoring the ... value
  ind<- match( names( cov.args), names(list(...) ) )
  cov.args = c(cov.args[is.na(ind)], list(...))
  # add lambda as a component to the starting parameters if missing
  if( is.null( cov.params.start$lambda)){
    cov.params.start <- c( cov.params.start, list( lambda= .5))
  }
  ########################################################################
  #  evaluate likelihood for a grid of theta on log scale
  # maximizing over lambda.
  #########################################################################
  # if range or starting value  for range is missing use quantiles of pairwise
  # distances among data.  
  cov.params.startTemp <-  cov.params.start 
  cov.params.startTemp$theta<- NULL
  
  if( is.null( gridTheta) ){
    if( is.null( cov.args$Distance)){
      pairwiseD<- dist(x)
    }
    else{
      pairwiseD<- do.call(cov.args$Distance, list(x))
      pairwiseD<- pairwiseD[col(pairwiseD) > row( pairwiseD) ]
    }
    theta.range<- quantile( pairwiseD, c(.02,.97))
    gridTheta  <- 10**seq( log10(theta.range[1]), log10(theta.range[2]), length.out=gridN )
  }
  # 
  par.grid<- list( theta= gridTheta )
  
  timeGrid<- system.time(
  MLEGrid<- mKrigMLEGrid(x, y,  
                         weights = weights, Z= Z, 
                      mKrig.args = mKrig.args,
                         cov.fun = cov.function, 
                       cov.args  = cov.args,
                        par.grid = par.grid, 
                           na.rm = na.rm,
                         verbose = verbose,
                            REML = REML,
                cov.params.start = cov.params.startTemp)
  )
  if( verbose){
    cat("mKrigMLEGrid summary", fill=TRUE)
     print(MLEGrid$summary)
  }
  #######################################################################
  # refine MLE for lambda and theta use the best value of theta from grid 
  # search if starting value not passed. 
  ########################################################################
  
    ind<- which.max( MLEGrid$summary[,"lnProfileLike.FULL"] )
    theta.start <-  par.grid$theta[ind]
    lambda.start<- MLEGrid$summary[ind,"lambda"] 
    cov.params.startTemp <- cov.params.start 
    # update starting values with results from grid search over theta
    cov.params.startTemp$lambda <- lambda.start
    cov.params.startTemp$theta  <- theta.start
    
  timeOptim<- system.time(
  MLEJoint <- mKrigMLEJoint(x, y, weights = weights, Z = Z,
                                            mKrig.args = mKrig.args,
                                               cov.fun = cov.function,
                                              cov.args = cov.args, 
                                      cov.params.start = cov.params.startTemp,
                                            optim.args = optim.args,
                                                abstol = abstol,
                                                 na.rm = na.rm,
                                               verbose = verbose,
                                                  REML = REML)
  )
  if( verbose){
    cat("mKrigMLEJoint summary", fill=TRUE)
    print( MLEJoint$summary)
  }
  ###########################################################################
  # evaluate likelihood on grid of log lambda with MLE for theta
  # NOTE lambda.profile = FALSE makes this work.
  ###########################################################################
  par.grid<- list( lambda = (10^(seq( -2,2,length.out=gridN)  ))*MLEJoint$pars.MLE["lambda"] )
  thetaMLE<- MLEJoint$pars.MLE["theta"]
  #par.grid<- list( theta= rep(thetaMLE, gridN) )
  cov.params.startTemp <- cov.params.start
  # do _not_ optimze over lambda (we are profiling)
  cov.params.startTemp$lambda <- NULL
  cov.params.startTemp$theta  <- thetaMLE
  timeProfile<- system.time(
  MLEProfileLambda <- mKrigMLEGrid(x, y,  weights = weights, Z= Z,
                                          cov.fun = cov.function, 
                                        cov.args  = cov.args,
                                       mKrig.args = mKrig.args,
                                         par.grid = par.grid, 
                                            na.rm = na.rm,
                                 cov.params.start = cov.params.startTemp,
                                             REML = REML,
                                          verbose = FALSE) 
  )
  timingTable<- cbind( timeGrid, timeOptim, timeProfile)
  return(
     list( summary= MLEJoint$summary, MLEGrid= MLEGrid, MLEJoint=MLEJoint, 
          MLEProfileLambda=MLEProfileLambda, call=match.call(), 
          timingTable= timingTable)
     )
}
