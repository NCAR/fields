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

MLESpatialProcess <- function(x, y, weights = rep(1, nrow(x)), Z = NULL,
                            mKrig.args = NULL,
                            cov.function = "stationary.cov", 
                            cov.args = list(Covariance = "Matern",
                                            smoothness = 1), 
                          lambda.guess = .5,
                           theta.start = NULL, 
                           theta.range = NULL,
                                 gridN = 20,
                            optim.args = NULL,
                                 na.rm = TRUE,
                               verbose = FALSE,
                               ...) {
  if( verbose){
    cat(" MLESpatialProcess extra arguments:" , full=TRUE)
     print( names( list(...)))
  }
  # combine  list(...) with cov.args and omit duplicates favoring the ... value
  ind<- match( names( cov.args), names(list(...) ) )
  cov.args = c(cov.args[is.na(ind)], list(...))
  #
  # if range or starting guess for range is missing use quantiles of pairwise
  # distances among data.  Use median pairwise distance as starting guess
  if( is.null( theta.range) ){
    if( is.null( cov.args$Distance)){
      pairwiseD<- dist(x)
    }
    else{
    pairwiseD<- do.call(cov.args$Distance, list(x))
    pairwiseD<- pairwiseD[col(pairwiseD) > row( pairwiseD) ]
    }
    theta.range<- quantile( pairwiseD, c(.02,.97))
    if (is.null(theta.start)) {
      theta.start<-  exp(mean( log(theta.range)))
    }
  }
  print( theta.range)
  #  evaluate likelihood for a grid of theta on log scale maximizing over lambda.
  # set all arguments for the optim function
  thetaGrid<- seq( theta.range[1], theta.range[2], length.out=gridN )
  par.grid<- list( theta= thetaGrid)
  MLEGrid<- mKrigMLEGrid(x, y,  weights = weights, Z= Z, 
                         mKrig.args = mKrig.args,
                         cov.fun = cov.function, 
                       cov.args  = cov.args,
                        par.grid = par.grid, 
                          lambda = lambda.guess, 
                  lambda.profile = TRUE, 
                           na.rm = na.rm,
                         verbose = verbose) 
  #refine MLE for lambda and theta
  if(is.null(optim.args)) {
    optim.args = list(method = "BFGS", 
                      control = list(fnscale = -1, parscale = c(0.5, 0.5), 
                                     ndeps = c(0.05,0.05)))
  }
  MLEJoint = mKrigMLEJoint(x, y, weights = weights,  Z = Z,
                           mKrig.args = mKrig.args,
                           cov.fun = cov.function,
                           cov.args = cov.args, 
                                          lambda.guess = lambda.guess, 
                                      cov.params.guess = list(theta=theta.start), 
                                            optim.args = optim.args,
                                                 na.rm = na.rm,
                                               verbose = verbose)
  # evaluate likelihood on grid of log lambda with MLE for theta
  #NOTE lambda.profile = FALSE makes this work.
  lambdaGrid<-   10^(seq( -2,1,,gridN))
  par.grid<- list( theta= rep(MLEJoint$pars.MLE[2], gridN) )
  MLEProfileLambda <- mKrigMLEGrid(x, y,  weights = weights, Z= Z,
                                          cov.fun = cov.function, 
                                        cov.args  = cov.args,
                                       mKrig.args = mKrig.args,
                                         par.grid = par.grid, 
                                           lambda = lambdaGrid, 
                                   lambda.profile = FALSE, 
                                            na.rm = na.rm,
                                          verbose = verbose) 
  return(
     list(MLEGrid= MLEGrid, MLEJoint=MLEJoint, 
          MLEProfileLambda=MLEProfileLambda, call=match.call() )
     )
}