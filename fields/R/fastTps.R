"fastTps" <- function(x, Y, m = NULL, p = NULL, aRange, 
                      lon.lat = FALSE, find.trA=FALSE,  REML=FALSE,
                      ...) {
  x <- as.matrix(x)
  d <- ncol(x)
  if (is.null(p)) {
    if (is.null(m)) {
      m <- max(c(2, ceiling(d/2 + 0.1)))
    }
    p <- (2 * m - d)
    if (p <= 0) {
      warning(" m is too small to satisfy thin plate spline crierion \n
                    you must have 2*m - dimension >0 \n
                    smoothness of Wendland set at k =2")
    }
  }
  # special arguments to send to the wendland covariance/taper function.
  # see nearest.dist for some explanation of 'method'
  method <- ifelse(!lon.lat, 
                   "euclidean", "greatcircle")
  cov.args <- list(        k = max(c(p,2)),
                           Dist.args = list(method=method),
                           aRange = aRange
  )
  
  object<-spatialProcess(x, Y, 
                         cov.function = "wendland.cov", 
                         mKrig.args = list( m = m), 
                         cov.args = cov.args,  REML=REML,
                         ... )
  
  object$call<- match.call()
  class(object) <- c( "fastTps", class(object))
  return( object)
}
