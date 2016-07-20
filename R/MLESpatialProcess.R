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

MLESpatialProcess <- function(x, y,
                              weights = rep(1, nrow(x)),
                          lambda.start = .5,
                           theta.start = NULL, 
                           theta.range = NULL,
                                 gridN = 20,
                          cov.function = "stationary.cov", 
                              cov.args = list(Covariance = "Matern", smoothness = 1), 
                                     Z = NULL,
                              Distance = "rdist",
                            optim.args = NULL,
                               verbose = FALSE,
                               ...) {
  if( verbose){
    cat("extra arguments:" , full=TRUE)
     print( names( list(...)))
  }
  # set the distance argument in cov.args if necessary
  if(supportsArg(cov.function, arg="Distance") && is.null(cov.args$Distance))
    cov.args = c(cov.args, list(Distance = Distance))
  #
  # if range or starting guess for range is missing use quantiles of pairwise
  # distances among data.  Use median pairwise distance as starting guess
  
  if( is.null( theta.range) ){
    pairwiseD<- do.call(Distance, list(x))
    theta.range<- quantile( pairwiseD[col(pairwiseD) > row( pairwiseD) ], c(.05,.95))
  }
  thetaGrid<- seq( theta.range[1], theta.range[2], length.out=gridN )
  par.grid<- list( theta= thetaGrid)
  
  obj<- mKrigCheckXY( x, y, weights, Z, na.rm=TRUE)
  
  #
  if (is.null(theta.start)) {
    pairwiseD<- do.call(Distance, list(x))
    theta.start <- median( pairwiseD[col(pairwiseD) > row( pairwiseD) ])
  }
#    
  
  # set all arguments for the optim function
  if(is.null(optim.args)) {
    optim.args = list(method = "BFGS", 
                      control = list(fnscale = -1, parscale = c(0.5, 0.5), 
                                     ndeps = c(0.05,0.05)))
  }
  # do parameter optimization
  
  MLEGrid<- mKrigMLEGrid(obj$x, obj$y,  weights = obj$weights, Z= obj$Z,
                         cov.fun = "stationary.cov", 
                         cov.args  = list(Covariance = "Matern", smoothness = 1.0),
                         par.grid = par.grid, 
                         lambda = lambda.start, 
                         lambda.profile = TRUE, 
                         verbose = verbose) 
  
  MLEJoint = do.call("mKrigMLEJoint", c(list(obj$x, obj$y,
                                               weights = obj$weights,
                                          lambda.guess = lambda.start, 
                                      cov.params.guess = list(theta=theta.start), 
                                               cov.fun = cov.function,
                                              cov.args = cov.args, 
                                            optim.args = optim.args,
                                                     Z = obj$Z,
                                               verbose = verbose),
                                                  list(...))
                                     )
  lambda.start<-   10^(seq( -2,1,,gridN))
  par.grid<- list( theta= rep(MLEJoint$MLEInfo$MLEJoint$pars.MLE[2], gridN) )
  MLEProfileLambda <- mKrigMLEGrid(obj$x, obj$y,  weights = obj$weights, Z= obj$Z,
                                          cov.fun = "stationary.cov", 
                                          cov.args  = list(Covariance = "Matern", smoothness = 1.0),
                                          par.grid = par.grid, 
                                          lambda = lambda.start, 
                                          lambda.profile = FALSE, 
                                          verbose = verbose) 
  return(
     list(MLEGrid= MLEGrid, MLEJoint=MLEJoint, 
          MLEProfileLambda=MLEProfileLambda, call=match.call() )
     )
}