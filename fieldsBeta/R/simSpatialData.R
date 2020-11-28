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
simSpatialData<- function(object,  M = 1, 
    verbose = FALSE) {
    # object is assummed to have the spatial model according to the
    # components in a spatialProcess or mKrig object
    # important variance parameters estimated from the data  
        tau2 <- (object$tau.MLE.FULL)^2
    #   set tau2 equal to zero if not there
        if( is.null(tau2)){ tau2 <- 0}
    #    
        sigma <- object$sigma.MLE.FULL
    #
    xUnique<- object$x
    if( any( duplicated(xUnique)) ){
      stop('Can not handle repeated
           prediction locations')}
    N.full <- nrow(xUnique)
    if (verbose) {
        cat("N.full", N.full, fill = TRUE)
    }
    if( N.full > 5000){
      cat("WARNING: Number of locations for  simulation is large ( >5000)
              this may take some time to compute or exhaust the memory.",
          fill=FALSE)
    }
    #
    # Sigma is full covariance at the data locations and at prediction points.
    # not to be confused with the lowercase tau that is the nugget variance
    # 
    Sigma <- sigma * do.call(object$cov.function.name, c(object$args, 
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
      h.data <- t(Schol) %*% matrix( rnorm(N.full*M), N.full,M)
    # value of simulated field at observations
      out<- h.data 
    # add measurement error (aka the "nugget")  
      if( tau2 > 0){
        nugget.error<- sqrt(tau2)*matrix( rnorm(N.full*M), N.full,M)/object$weights
        out<- out + nugget.error
      }
    #
      return( out)
}
