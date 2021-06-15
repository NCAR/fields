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

mKrigMLEGrid <- function(x, y, weights = rep(1, nrow(x)), Z = NULL,
                       mKrig.args = NULL,
                          cov.function = "stationary.cov", 
                         cov.args = NULL,
                           na.rm = TRUE, 
                         par.grid = NULL, 
                           reltol = 1e-06,
                             REML = FALSE,
                             GCV  = FALSE,
                       optim.args = NULL,
                 cov.params.start = NULL,
                          verbose = FALSE) {
  if( na.rm){
    obj<- mKrigCheckXY(x, y, weights, Z, na.rm)
    x<- obj$x
    y<- obj$y
    weights<- obj$weights
    Z<- obj$Z
  }
  #check which optimization options the covariance function supports
  #precompute distance matrix if possible so it only needs to be computed once
  supportsDistMat = supportsArg(cov.function, "distMat")
  #precompute distance matrix if possible so it only needs to be computed once
  if(supportsDistMat) {
    #Get distance function and arguments if available
    #If user left all distance settings NULL, use rdist with compact option.
    #Use rdist function by default in general.
    #
   
    if(is.null(cov.args$Distance)) {
      cov.args$Distance  <-  "rdist"
      cov.args$Dist.args <- list(compact=TRUE)
    }
    cov.args$distMat<-do.call(cov.args$Distance, c( list(x), cov.args$Dist.args) )
    cov.args$onlyUpper<- TRUE
    }
# find NG --  number of parameters to try
  if( is.null(par.grid) ){ 
    stop("no par.grid ")
  }
  par.grid <- data.frame(par.grid)
  NG <- nrow(par.grid)
  # output object to summarize results
  summary <-NULL
  # begin loop over covariance parameters and either evaluate at a grid of  lambda values
  # or use these as start values for optimization.
  for (k in 1:NG) {
    cov.args.temp <- as.list( par.grid[k, ])
    names(cov.args.temp) <- names( par.grid)
    currentCov.args<- c( cov.args.temp, cov.args) 
    if( verbose){
      
      cat( "********grid value: " , k, fill=TRUE)
      cat( "Cov args", names(currentCov.args ), fill=TRUE, sep=", ")
      cat("value  aRange", fill=TRUE)
      print( currentCov.args$aRange)
      cat("value  lambda", fill=TRUE)
      print( currentCov.args$lambda)
      
      cat("start values", fill=TRUE)
      print( cov.params.start)
      cat("optim.args", fill=TRUE)
      print( optim.args)
    }
  
    MLEfit0 <- mKrigMLEJoint(x, y, 
                                  weights = weights, Z=Z, 
                                  cov.function = cov.function,
                               optim.args = optim.args,
                                 cov.args = currentCov.args,
                                    na.rm = na.rm,
                               mKrig.args = mKrig.args,
                                     REML = REML,
                                     GCV  = GCV,
                                   reltol = reltol,
                         cov.params.start = cov.params.start,
                                  verbose = verbose)
     summary <- rbind( summary, MLEfit0$summary)
  }
  summary<- cbind( summary, par.grid)
  if( REML){
    indMax<- which.max( summary[,"lnProfileLike.FULL"]) 
  } 
  else{
    indMax<- which.max( summary[,"lnProfileREML.FULL"]) 
  }
  if( verbose){
    cat("summary from mKrigMLEGrid:", fill=TRUE)
    print(summary )
  }
  return(list(summary = summary, par.grid = par.grid,
              call = match.call(), indMax= indMax )
         )
}
