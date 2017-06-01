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

mKrigMLEGrid <- function(x, y, weights = rep(1, nrow(x)), Z = NULL,
                       mKrig.args = NULL,
                          cov.fun = "stationary.cov", 
                         cov.args = NULL,
                           na.rm = TRUE, 
                         par.grid = NULL, 
                           lambda = NULL, 
                   lambda.profile = TRUE, 
               relative.tolerance = 1e-04,
                             REML = FALSE,
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
  supportsDistMat = supportsArg(cov.fun, "distMat")
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

  lnProfileLike.max <- -1e+20
# find NG --  number of parameters to try
  par.grid <- data.frame(par.grid)
  if (nrow(par.grid) == 0) {
    NG<- ifelse(is.null(lambda), 1, length( lambda)) 
  }
  else {
    NG <- nrow(par.grid)
  }
  lambda.best <- NA
  # default for lambda is 1.0 for first value and exp(llambda.opt) for subsequent ones
  # this is controlled by NAs for lambda starting values.
  if (is.null(lambda)) {
    lambda <- rep(NA, NG)
  }
  # output matrix to summarize results
  summary <- matrix(NA, nrow = NG, ncol = 8)
  
  # default starting value for lambda is .5 or log lambda is 0
  lambda.opt <- .5
  optim.counts <- c(NA, NA)
  lnLike.eval <- list()
  # Define the objective function as a tricksy call to mKrig
  # if Y is a matrix of replicated data sets use the log likelihood for the complete data sets
  #
  # begin loop over covariance arguments
  lnLike.eval<- list()
  for (k in 1:NG) {
    lambda.start <- ifelse(is.na(lambda[k]), lambda.opt, (lambda[k]))
    # list of covariance arguments from par.grid with right names (some R arcania!)
    # note that this only works because 1) temp.fn will search in this frame for this object
    # par.grid has been coerced to a data frame so one has a concept of a row subscript.
    cov.args.temp <- as.list(par.grid[k, ])
    names(cov.args.temp) <- names(par.grid)
    currentCov.args<- c(cov.args.temp, cov.args) 
    # optimize over lambda if lambda.profile is TRUE
    optim.args = list(method = "BFGS", 
                      control = list(fnscale = -1, parscale = c(0.5), 
                                     ndeps = c(0.05)))
    if (lambda.profile) {
      # set up matrix to store evaluations from within optim
    MLEfit0 <- mKrigMLEJoint(x, y, weights=weights, Z=Z, 
                             lambda.start = lambda.start, 
                         cov.params.start = NULL, 
                                  cov.fun = cov.fun,
                               optim.args = optim.args,
                                 cov.args = currentCov.args,
                                    na.rm = na.rm,
                               mKrig.args = mKrig.args,
                                     REML = REML,
                                  verbose = verbose)
    lnLike.eval<- c( lnLike.eval, list(MLEfit0$lnLike.eval))
    lambda.opt<- MLEfit0$pars.MLE[1]
    }
    else {
      # no refinement for lambda so just save the the 'start' value as final one.
      lambda.opt <- lambda.start
    }
    
# final fit at optimal value 
#    (or starting value if not refinement/maximization for lambda)
    obj <- do.call("mKrig", c(
      list(x = x, y = y, weights = weights, Z = Z, na.rm = na.rm),
                                    mKrig.args,
       list(lambda=lambda.opt),
       list( cov.fun= cov.fun, cov.args = currentCov.args)
      )
      )
    nameCriterion<- ifelse( !REML,
                            "lnProfileLike.FULL",
                            "lnProfileREML.FULL" )
    if (obj[[nameCriterion]] > lnProfileLike.max) {
      lnProfileLike.max <- obj$lnProfileLike.FULL
      cov.args.MLE <- cov.args.temp
      lambda.best <- lambda.opt
    }
    
# save results of the kth covariance model evaluation
    summary[k, 1:8] <- c(obj$eff.df, obj[[nameCriterion]], 
                         obj$GCV, obj$sigma.MLE.FULL, obj$rho.MLE.FULL, lambda.opt, 
                         optim.counts)
    dimnames(summary) <- list(NULL, c("EffDf",nameCriterion , 
                                      "GCV", "sigma.MLE", "rho.MLE", "lambda.MLE", "counts eval", 
                                      "counts grad"))
    if (verbose) {
      cat("Summary: ", k, summary[k, 1:8], fill = TRUE)
    }
  }
  return(list(summary = summary, par.grid = par.grid, cov.args.MLE = cov.args.MLE, 
              lambda.best = lambda.best, lambda.MLE = lambda.best, 
              call = match.call(), lnLike.eval = lnLike.eval)
         )
}
