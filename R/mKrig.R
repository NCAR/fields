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
mKrig <- function(x, y, weights=rep(1, nrow(x)), Z = NULL,
                  cov.function="stationary.cov", 
                  cov.args = NULL, lambda = 0, m = 2, 
                  chol.args = NULL, find.trA = TRUE, NtrA = 20, 
                  iseed = 123, llambda = NULL, na.rm=FALSE, ...) {
  
  # pull extra covariance arguments from ...  and overwrite
  # any arguments already named in cov.args
  ind<- match( names( cov.args), names(list(...) ) )
  cov.args = c(cov.args[is.na(ind)], list(...))
  #
  #If cov.args$find.trA is true, set onlyUpper to FALSE (onlyUpper doesn't
  #play nice with predict.mKrig, called by mKrig.trace)
  #
  if(find.trA == TRUE && supportsArg(cov.function, "onlyUpper"))
    cov.args$onlyUpper= FALSE
  if(find.trA == TRUE && supportsArg(cov.function, "distMat"))
    cov.args$distMat= NA
  
  if (!is.null(llambda)) {
    lambda <- exp(llambda)
  }
  # see comments in Krig.engine.fixed for algorithmic commentary
  #
  # check for duplicate x's.
  # stop if there are any
  if (any(duplicated(cat.matrix(x)))) {
    stop("locations are not unique see help(mKrig) ")
  }
  # next function also omits NAs from x,y,weights, and Z  if na.rm=TRUE.
  object<- mKrigCheckXY( x, y, weights, Z, na.rm = na.rm)
  # create fixed part of model as m-1 order polynomial
  # NOTE: if m==0 then fields.mkpoly returns a NULL to 
  # indicate no polynomial part.
  Tmatrix <- cbind(fields.mkpoly(object$x, m), object$Z)
  # set some dimensions
    np <- nrow(object$x)
  if( is.null(Tmatrix) ){
    nt<- 0
    }
  else{
    nt<- ncol(Tmatrix) 
  }
  if( is.null(object$Z)){
    nZ<- 0
  }
  else{
    nZ<- ncol(object$Z)
  }
  ind.drift <- c(rep(TRUE, (nt - nZ)), rep(FALSE, nZ)) 
  
  # as a place holder for reduced rank Kriging, distinguish between
  # observations locations and  the locations to evaluate covariance.
  # (this is will also allow predict.mKrig to handle a Krig object)
  object$knots <- object$x
  # covariance matrix at observation locations
  # NOTE: if cov.function is a sparse constuct then Mc will be sparse.
  # see e.g. wendland.cov
  Mc <- do.call(cov.function, c(cov.args, list(x1 = object$knots, x2 = object$knots)))
  #
  # decide how to handle the pivoting.
  # one wants to do pivoting if the matrix is sparse.
  # if Mc is not a matrix assume that it is in sparse format.
  #
  sparse.flag <- !is.matrix(Mc)
  #
  # set arguments that are passed to cholesky
  #
  if (is.null(chol.args)) {
    chol.args <- list(pivot = sparse.flag)
  }
  else {
    chol.args <- chol.args
  }
  # quantify sparsity of Mc for the mKrig object
  nzero <- ifelse(sparse.flag, length(Mc@entries), np^2)
  # add diagonal matrix that is the observation error Variance
  # NOTE: diag must be a overloaded function to handle sparse format.
  if (lambda != 0) {
    if(! sparse.flag)
      invisible(.Call("addToDiagC", Mc, as.double(lambda/object$weights), nrow(Mc)))
    else
      diag(Mc) = diag(Mc) + lambda/object$weights
  }
  #  MARK LINE Mc
  # At this point Mc is proportional to the covariance matrix of the
  # observation vector, y.
  #
  # cholesky decoposition of Mc
  # do.call used to supply other arguments to the function
  # especially for sparse applications.
  # If chol.args is NULL then this is the same as
  #              Mc<-chol(Mc), chol.args))
  Mc <- do.call("chol", c(list(x = Mc), chol.args))
 
  lnDetCov <- 2 * sum(log(diag(Mc)))
  
  #
  # start linear algebra to find estimates and likelihood
  # Note that all these expressions make sense if y is a matrix
  # of several data sets and one is solving for the coefficients
  # of all of these at once. In this case d.coef and c.coef are matrices
  #
  if( !is.null(Tmatrix)){
  # Efficent way to multply inverse of Mc times the Tmatrix  
  VT <- forwardsolve(Mc, x = Tmatrix, k=ncol(Mc), transpose = TRUE, upper.tri = TRUE)
  qr.VT <- qr(VT)
  
  # now do generalized least squares for d
    d.coef <- as.matrix(qr.coef(qr.VT, forwardsolve(Mc, transpose = TRUE, 
                                                  object$y, upper.tri = TRUE)))
    resid<-  object$y - Tmatrix %*% d.coef
  # GLS covariance matrix for fixed part.
    Rinv <- solve(qr.R(qr.VT))
    Omega <- Rinv %*% t(Rinv)
#    
#  Omega is  solve(t(Tmatrix)%*%solve( Sigma)%*%Tmatrix)
# proportional to fixed effects covariance matrix.    
#  Sigma = cov.function( x,x) + lambda/object$weights
#  this is proportional to the covariance matrix for the GLS estimates of
#  the fixed linear part of the model. 
#     
    R2diag<-  diag( qr.R(qr.VT) )^2
    lnDetOmega<- -1* sum( log(R2diag) ) 
  }
  else{
# much is set to NULL because no fixed part of model    
    nt<- 0
    resid<- object$y
    Rinv<- NULL
    Omega<- NULL
    qr.VT<- NULL
    d.coef<- NULL
    lnDetOmega <- 0
  }
  # and now find c.
  #  the coefficents for the spatial part.
  # if linear fixed part included resid as the residuals from the 
  # GLS regression.
  c.coef <- as.matrix(forwardsolve(Mc, transpose = TRUE,
                                   resid, upper.tri = TRUE))
  # save intermediate result this is   t(y- T d.coef)( M^{-1}) ( y- T d.coef)
  quad.form <- c(colSums(as.matrix(c.coef^2)))
  # find c coefficients
  c.coef <- as.matrix(backsolve(Mc, c.coef))
  # MLE estimate of rho and sigma
  #    rhohat <- c(colSums(as.matrix(c.coef * y)))/(np - nt)
  # NOTE if y is a matrix then each of these are vectors of parameters.
  rho.MLE <- quad.form/np
  rhohat <- c(colSums(as.matrix(c.coef * object$y)))/np
  shat.MLE <- sigma.MLE <- sqrt(lambda * rho.MLE)
  # the  log profile likehood with  rhohat  and  dhat substituted
  # leaving a profile for just lambda.
  # NOTE if y is a matrix then this is a vector of log profile
  # likelihood values.
  lnProfileLike <- (-np/2 - log(2 * pi) * (np/2) - (np/2) * 
                      log(rho.MLE) - (1/2) * lnDetCov)
  # see section 4.2 handbook of spatial statistics (Zimmermanchapter)
  lnProfileREML <-  lnProfileLike + (1/2) * lnDetOmega
  rho.MLE.FULL <- mean(rho.MLE)
  sigma.MLE.FULL <- sqrt(lambda * rho.MLE.FULL)
  # if y is a matrix then compute the combined likelihood
  # under the assumption that the columns of y are replicated
  # fields
  lnProfileLike.FULL <- sum((-np/2 - log(2 * pi) * (np/2) - 
                               (np/2) * log(rho.MLE.FULL) 
                               - (1/2) * lnDetCov)
                            )
  lnProfileREML.FULL <- sum((-np/2 - log(2 * pi) * (np/2) - 
                               (np/2) * log(rho.MLE.FULL) 
                             - (1/2) * lnDetCov
                             + (1/2) * lnDetOmega  )
                            )
  
  #
  # return coefficients and   include lambda as a check because
  # results are meaningless for other values of lambda
  # returned list is an 'object' of class mKrig (micro Krig)
  # also save the matrix decompositions so coefficients can be
  # recalculated for new y values.  Make sure onlyUpper and 
  # distMat are unset for compatibility with mKrig S3 functions
  if(!is.null(cov.args$onlyUpper))
    cov.args$onlyUpper = FALSE
  if(!is.null(cov.args$distMat))
    cov.args$distMat = NA
  object <- c( object, list( 
               d = d.coef, c = c.coef, nt = nt, np = np, 
              lambda.fixed = lambda, 
              cov.function.name = cov.function, 
              args = cov.args, m = m, chol.args = chol.args, call = match.call(), 
              nonzero.entries = nzero, shat.MLE = sigma.MLE, sigma.MLE = sigma.MLE, 
              rho.MLE = rho.MLE, rhohat = rho.MLE, lnProfileLike = lnProfileLike, 
              rho.MLE.FULL = rho.MLE.FULL, sigma.MLE.FULL = sigma.MLE.FULL, 
              lnProfileLike.FULL = lnProfileLike.FULL,
              lnProfileREML.FULL =  lnProfileREML.FULL,
              lnProfileREML =  lnProfileREML,
              lnDetCov = lnDetCov, lnDetOmega = lnDetOmega,
              quad.form = quad.form, Omega = Omega,lnDetOmega=lnDetOmega,
              qr.VT = qr.VT, 
              Mc = Mc, Tmatrix = Tmatrix, ind.drift = ind.drift, nZ = nZ)
  )
  #
  # find the residuals directly from solution
  # to avoid a call to predict
  object$residuals <- lambda * c.coef/object$weights
  object$fitted.values <- object$y - object$residuals
  # estimate effective degrees of freedom using Monte Carlo trace method.
  if (find.trA) {
    object2 <- mKrig.trace(object, iseed, NtrA)
    object$eff.df <- object2$eff.df
    object$trA.info <- object2$trA.info
    object$GCV <- (sum(object$residuals^2)/np)/(1 - object2$eff.df/np)^2
    if (NtrA < np) {
      object$GCV.info <- (sum(object$residuals^2)/np)/(1 - object2$trA.info/np)^2
    }
    else {
      object$GCV.info <- NA
    }
  }
  else {
    object$eff.df <- NA
    object$trA.info <- NA
    object$GCV <- NA
  }
  class(object) <- "mKrig"
  return(object)
}
