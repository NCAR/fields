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
#

predict.fastTps <- function(object, xnew = NULL, grid.list=NULL,
                            ynew = NULL,  derivative = 0,
                            Z = NULL, drop.Z = FALSE, just.fixed = FALSE, xy=c(1,2), ...)
{
  # the main reason to pass new args to the covariance is to increase
  # the temp space size for sparse multiplications
  # other optional arguments from mKrig are passed along in the
  # list object$args
  cov.args <- list(...)
  # predict using grid.list or as default observation locations
  if( !is.null(grid.list)){
    xnew<- make.surface.grid( grid.list)
  }
  if( is.null(xnew) ) {
    xnew <- object$x
  }   
  if (!is.null(ynew)) {
    coef.hold <- mKrig.coef(object, ynew)
    c.coef <- coef.hold$c.coef
    beta <- coef.hold$beta
  }
  else {
    c.coef <- object$c.coef
    beta <- object$beta
  }
  # fixed part of the model this a polynomial of degree m-1
  # Tmatrix <- fields.mkpoly(xnew, m=object$m)
  #
  if (derivative == 0){
    if (drop.Z | object$nZ == 0) {
      # just evaluate polynomial and not the Z covariate
      temp1 <- fields.mkpoly(xnew, m = object$m) %*% beta[object$ind.drift, ]
    }
    else{
      if( is.null(Z)) {
        Z <- object$Tmatrix[, !object$ind.drift]
      }
      temp1 <- cbind(fields.mkpoly(xnew, m = object$m), Z) %*% beta
    }
  }
  else{
    if (!drop.Z & object$nZ > 0) {
      stop("derivative not supported with Z covariate included")
    }
    temp1 <- fields.derivative.poly(xnew, m = object$m, beta[object$ind.drift, 
    ])
  }
  if (just.fixed) {
    return(temp1)
  }
  
  useFORTRAN<-  (ncol(object$x)==2) & (object$args$k == 2) & (derivative==0) & (!is.null(grid.list))
  
  
  # add nonparametric part. 
  # call FORTRAN under a specific case  
  if( useFORTRAN){
    
    temp2<- multWendlandGrid(grid.list, object$knots, delta=object$args$aRange, c.coef, xy=xy)
    
  }
  else{
    temp2 <- do.call(object$cov.function.name, c(object$args, 
                                                 list(x1 = xnew, x2 = object$knots, C = c.coef, derivative = derivative), 
                                                 cov.args))
  }  
  # add two parts together
  return((temp1 + temp2))
}
