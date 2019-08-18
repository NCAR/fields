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
sim.spatialProcess<- function(object, xp,  M = 1, 
    verbose = FALSE, ...) {
        sigma2 <- (object$sigma.MLE.FULL)^2
        rho <- object$rho.MLE.FULL
        xp<- as.matrix( xp)
    #
    # check for unique rows of  data locations
    if( any(duplicated(object$x)) ){
        stop("data locations should be unique")
    }
    #
    # set up various sizes of arrays
    m <- nrow(xp)
    n<- nrow( object$x)
    N <- length(object$y)
    if (verbose) {
        cat("m,n,N", m,n, N, fill = TRUE)
    }
    #
    # find indices of all rows of xp that correspond to rows of
    # xM and then collapse x to unique rows.
    if( any( duplicated(object$x)) ){
      stop('Can not handle replicated locations')}
    if( any( duplicated(xp)) ){
      stop('Can not handle repeated
           prediction locations')}
   # 
    x<- as.matrix(rbind( object$x, xp))
    rep.x.info <- fields.duplicated.matrix(x)
    # find uniuqe locations. 
    ind<- !duplicated(rep.x.info)
    xUnique <- as.matrix(x[ind, ])
    if (verbose) {
        cat('full x and predicted locations without duplicates ', fill = TRUE)
        print(xUnique)
    }
    N.full <- nrow(xUnique)
    if (verbose) {
        cat("N.full", N.full, fill = TRUE)
    }
    if( N.full > 5000){
      cat("WARNING: Number of locations for conditional simulation is large ( >5000)
              this may take some time to compute or exhaust the memory.",
          fill=FALSE)
    }
    # these give locations in x matrix to reconstruct xp matrix
    xp.ind <- rep.x.info[(1:m) + n]
    if (verbose) {
        print(N.full)
        print(xUnique)
    }
    if (verbose) {
        cat("reconstruction of xp from collapsed locations", 
            fill = TRUE)
        print(xUnique[xp.ind, ])
    }
    #
    # Sigma is full covariance at the data locations and at prediction points.
    #
    
    Sigma <- rho * do.call(object$cov.function.name, c(object$args, 
        list(x1 = xUnique, x2 = xUnique)))
    #
    # square root of Sigma for simulating field
    # Cholesky is fast but not very stable.
    #
    # the following code line is similar to chol(Sigma)-> Scol
    # but adds possible additional arguments controlling the Cholesky
    # from the mKrig object.
    # x has has been winnowed down to unique rows so that 
    # Sigma has full rank. 
    #
    Schol <- do.call("chol", c(list(x = Sigma), object$chol.args))
    #
    # output matrix to hold results
    out <- matrix(NA, ncol = m, nrow = M)
    #
    # find conditional mean field from initial fit
    # (these are added at the predict step).
    #
     h.hat <- predict(object, xnew=xp,  ...)
    #
    # NOTE: fixed part of model (null space) need not be simulated
    # because the  estimator is unbiased for this part.
    # the variability is still captured because the fixed part
    # is still estimated as part of the predict step below
    # create synthetic data
    for (k in 1:M) {
        # simulate full field
        h <- t(Schol) %*% rnorm(N.full)
        # value of simulated field at observations
        h.data <- h[1:n]
        #
        y.synthetic <- h.data + sqrt(sigma2/object$weights)*rnorm(n)
        # predict at xp using these data
        # and subtract from 'true' value, 
        # note that true values of field have to be expanded in the
        # case of common locations between object$x and xp.
        h.true <- (h[xp.ind])
#        temp.error <- predict(object, xnew=xp, ynew = y.synthetic, 
#                              Z=Zp, ...) - h.true
        temp.error <- predict(object, xnew=xp, ynew = y.synthetic, 
                                         ...) - h.true
        # add the error to the actual estimate  (conditional mean)
        out[k, ] <- h.hat + temp.error 
    }
    out
}
