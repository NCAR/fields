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
mKrigCheckXY <- function(x, y,  weights, Z, na.rm) 
    {
  #
  # check for missing values in y or X.
  #
  # save logical indicating where there are NA's
  # and check for NA's
  #
  ind <- is.na(y)
  if (any(ind) & !na.rm) {
    stop("Need to remove missing values or use: na.rm=TRUE in the call")
  }
  #
  # coerce x to be a matrix
  x <- as.matrix(x)
  #
  # coerce y to be a matrix
  #
  y <- as.matrix(y)
  #
  #default weights ( reciprocal variance of errors).
  #
  if (is.null(weights)) 
    weights <- rep(1, nrow(y))
  #
  # check that dimensions agree
  #
  if (nrow(y) != nrow(x)) {
    stop(" length of y and number of rows of x differ")
  }
  if (nrow(y) != length(weights)) {
    stop(" length of y and weights differ")
  }
  #  if Z is not NULL coerce to be  a matrix
  # and check  # of rows
  if (!is.null(Z)) {
    if (!is.matrix(Z)) {
      Z <- as.matrix(Z)
    }
    if (nrow(y) != nrow(Z)) {
      stop(" number of rows of y and number of rows of Z differ")
    }
  }
  # if NAs can be removed then remove them and warn the user
  if (na.rm) {
    ind <- is.na(y)
    if(all(ind)){
      stop("Oops! All y values are missing!")
    }
    if (any(ind)) {
      y <- y[!ind]
      x <- as.matrix(x[!ind, ])
      if (!is.null(Z)) {
        Z <- as.matrix(Z[!ind, ])
      }
      weights <- weights[!ind]
    }
  }
  #
  # check for NA's in x matrix -- there should not be any !
  if (any(c(is.na(x)))) {
    stop(" NA's in x matrix")
  }
  #
  # check for NA's in Z matrix
  if (!is.null(Z)) {
    if (any(c(is.na(Z)))) {
      stop(" NA's in Z matrix")
    }
  }
 
  # save x, weights  and y w/o NAs
  N <- length(y)
  return(list(N = N, y = y, x = x, weights = weights, Z = Z, 
              NA.ind = ind, nreps = ncol( y) )
         )
}
