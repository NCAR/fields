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

mKrig.trace <- function(object, iseed, NtrA) {
# do not reset the random seed if NA. 
  if( !is.na(iseed)){
    set.seed(iseed)
  }
    # if more MonteCarlo samples > number of data points just
    # find A exactly using  np  calls to predict.
    np<- object$np
    if (NtrA >= object$np) {
        Ey <- diag(1, np)
        NtrA <- np
        hold <- diag(predict.mKrig(object, ynew = Ey, collapseFixedEffect=FALSE))
        trA.info<- NA
        trA.est <- sum(hold)
    }
    else {
        # if fewer tests then use random trace method
        # find fitted.values  for iid N(0,1) 'data' to calculate the
        # the Monte Carlo estimate of tr A(lambda)
        # basically repeat the steps above but take some
        # short cuts because we only need fitted.values
        # create random normal 'data'
        Ey <- matrix(rnorm(np * NtrA), nrow = np, 
            ncol = NtrA)
        trA.info <- colSums(Ey * (predict.mKrig(object, ynew = Ey,
                                    collapseFixedEffect=FALSE)))
        trA.est <- mean(trA.info)
    }
    if (NtrA < np) {
     MSE<-(sum(object$residuals^2)/np) 
     GCV <-       MSE/(1 - trA.est /np)^2
     GCV.info <- MSE/( 1 - trA.info/np)^2
    }
    else{
    	GCV<- NA
    	GCV.info <- NA
    }	
    return(
    list(trA.info = trA.info, eff.df = trA.est,
             GCV= GCV, GCV.info=GCV.info)
    )
}

mKrig.coef <- function(object, y, collapseFixedEffect=TRUE) {
    # given new data y and the matrix decompositions in the
    # mKrig object find coficients beta  and c.
    # beta are the coefficients for the fixed part
    # in this case hard coded for a low order polynomial
    # c are coefficients for the basis functions derived from the
    # covariance function.
    #
    # see mKrig itself for more comments on the linear algebra
    #
    # Note that all these expressions make sense if y is a matrix
    # of several data sets and one is solving for the coefficients
    # of all of these at once. In this case beta and c.coef are matrices
    #
    # generalized least squares for d
    if( any(is.na(y))){
    	stop("mKrig can not omit missing values in observation vecotor")
    }
   if( object$nt>0){
    beta <- as.matrix(qr.coef(object$qr.VT, forwardsolve(object$Mc, 
        transpose = TRUE, y, upper.tri = TRUE)))
    betaMeans<- rowMeans( beta)
    if( collapseFixedEffect){ 
      beta<- matrix( betaMeans, ncol=ncol(beta), nrow= nrow( beta))
    }
    #  residuals from subtracting off fixed part
    #  of model as m-1 order polynomial
    resid <- y - object$Tmatrix %*% beta
   }
  else{
    beta<- NULL
    resid <- y
  }
    # and now find c.
    c.coef <- forwardsolve(object$Mc, transpose = TRUE, resid, 
        upper.tri = TRUE)
    c.coef <- as.matrix(backsolve(object$Mc, c.coef))
    out <- list( beta = (beta), c.coef = c.coef )
    return(out)
}






