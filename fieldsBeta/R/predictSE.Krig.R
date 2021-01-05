
"predictSE.Krig" <- function(object, x = NULL, cov = FALSE, 
                             verbose = FALSE, ...) {
  #
  # name of covariance function
  call.name <- object$cov.function.name
  #
  # default is to predict at data x's
  if (is.null(x)) {
    x <- object$x
  }
  x <- as.matrix(x)
  if (verbose) {
    print(x)
  }
  xraw <- x
  # transformations of x values used in Krig
  # NOTE knots are already scaled in Krig object
  xc <- object$transform$x.center
  xs <- object$transform$x.scale
  x <- scale(x, xc, xs)
  #
  # scaled unique observation locations.
  xM <- object$xM
  # find marginal variance before transforming x.
  
  temp.sd <- 1
  
  if (verbose) {
    print(temp.sd)
  }
  # Default is to use parameters in best.model
  lambda <- object$best.model[1]
  sigma2 <- object$best.model[3]
  tau2 <- object$best.model[2]
  nx <- nrow(xM)
  wght.vec <- t(Krig.Amatrix(object, xraw, lambda, eval.correlation.model = FALSE, 
                             ...))
  if (verbose) {
    cat("wght.vector", fill = TRUE)
    print(wght.vec)
  }
  #var( f0 - yhat)=    var( f0) -  cov( f0,yhat) - cov( yhat, f0) +  cov( yhat)
  #               =      temp0  - temp1 - t( temp1)       + temp2
  #
  # if off diagonal weight matrix is passed then
  # find  inverse covariance matrix
  # otherwise just create this quickly from diagonal weights
  #
  Wi <- Krig.make.Wi(object)$Wi
  # find covariance of data
  if (object$nondiag.W) {
    Cov.y <- sigma2 * do.call(call.name, c(object$args, list(x1 = xM, 
                                                             x2 = xM))) + tau2 * Wi
  }
  else {
    #     this is one case where keeping diagonal
    #     matrix as a vector will not work.
    Cov.y <- sigma2 * do.call(call.name, c(object$args, list(x1 = xM, 
                                                             x2 = xM))) + tau2 * diag(Wi)
  }
  if (!cov) {
    # find diagonal elements of covariance matrix
    # now find the three terms.
    # note the use of an element by element multiply to only get the
    # diagonal elements of the full
    #  prediction covariance matrix.
    #
    temp1 <- sigma2 * colSums(wght.vec * do.call(call.name, 
                                                 c(object$args, list(x1 = xM, x2 = x))))
    temp2 <- colSums(wght.vec * (Cov.y %*% wght.vec))
    #
    # find marginal variances -- trival in the stationary case!
    # Note that for the case of the general covariances
    # as radial basis functions (RBFs) temp0 should be zero.
    # Positivity results from the generalized divided difference
    # properties of RBFs.
    temp0 <- sigma2 * do.call(call.name, c(object$args, list(x1 = x, 
                                                             marginal = TRUE)))
    #
    temp <- temp0 - 2 * temp1 + temp2
    #
    return(sqrt(temp * temp.sd^2))
  }
  else {
    #
    # find full covariance matrix
    #
    temp1 <- sigma2 * t(wght.vec) %*% do.call(call.name, c(object$args, 
                                                           list(x1 = xM, x2 = x)))
    #
    temp2 <- t(wght.vec) %*% Cov.y %*% wght.vec
    #
    temp0 <- sigma2 * do.call(call.name, c(object$args, list(x1 = x, 
                                                             x2 = x)))
    #
    temp <- temp0 - t(temp1) - temp1 + temp2
    temp <- t(t(temp) * temp.sd) * temp.sd
    #
    return(temp)
  }
}
