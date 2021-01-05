"predictSE.mKrig" <- function(object, xnew = NULL, 
                              Z = NULL, verbose = FALSE, drop.Z = FALSE, ...) {
  #
  # name of covariance function
  call.name <- object$cov.function.name
  #
  # default is to predict at data x's
  if (is.null(xnew)) {
    xnew <- object$x
  }
  if ((!drop.Z) & !is.null(object$Z)) {
    Z <- object$Z
  }
  xnew <- as.matrix(xnew)
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
  }
  if (verbose) {
    print(xnew)
    print(Z)
  }
  objectSummary<- object$summary
  lambda <- objectSummary["lambda"]
  sigma2 <- objectSummary["sigma2"]
  tau2 <- objectSummary["tau"]^2
  if (verbose) {
    print(c(lambda, sigma2, tau2))
  }
  k0 <- do.call(call.name, c(object$args, list(x1 = object$x, 
                                               x2 = xnew)))
  # fixed effects matrox includes both spatial drift and covariates.
  if (!drop.Z) {
    t0 <- t(cbind(fields.mkpoly(xnew, m = object$m), Z))
  }
  else {
    stop(" drop.Z not supported")
  }
  #
  # old form based on the predict function
  #   temp1 <-  sigma2*(t0%*% object$Omega %*%t(t0)) -
  #          sigma2*predict( object, y= k0, x=x) -
  #          sigma2*predict( object, y= k0, x=x, just.fixed=TRUE)
  
  # alternative formula using the beta and c.coef coefficients directly.
  # collapseFixedEffect=FALSE because
  # we want the "fixed effect" computation
  # to be done separately for each column of k0
  hold <- mKrig.coef(object, y = k0, collapseFixedEffect=FALSE)
  temp1 <- sigma2 * (colSums(t0 * (object$Omega %*% t0)) - 
                       colSums((k0) * 
                       hold$c.coef) - 2 * colSums(t0 * hold$beta))
  # find marginal variances -- trival in the stationary case!
  temp0 <- sigma2 * do.call(call.name,
                            c(object$args, list(x1 = xnew, marginal = TRUE)))
  # Add marginal variance to part from estimate
  temp <- temp0 + temp1
  # return square root as the standard error in units of observations.
  return(sqrt(temp))
}
