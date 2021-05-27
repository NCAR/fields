# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2021
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

"predictSurfaceSE"<- function( object,...){
  UseMethod("predictSurfaceSE")
}

"predictSurfaceSE.default" <- 
  function(object, grid.list = NULL, 
           extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
           xy = c(1,2),  verbose = FALSE,
           ZGrid=NULL, just.fixed=FALSE,  ...) {
    
    # create a default grid if it is not passed    
    if (is.null(grid.list)) {
      # NOTE: 
      # without grid.list
      # default is 80X80 grid on first two variables
      # rest are set to median value of the x's
      grid.list <- fields.x.to.grid(object$x, nx = nx, ny = ny, 
                                    xy = xy)
    }
    # do some checks on Zgrid if passed and also reshape as a matrix
    # rows index grid locations and columns  are the covariates
    # (as Z in predict).
    # if ZGrid is NULL just returns that back
    
    Z<- unrollZGrid( grid.list, ZGrid) 
    
    # Convert grid.list to the full set of locations
    xg <- make.surface.grid(grid.list)
    # NOTE: the predict function called will need to do some internal  the checks
    # whether the evaluation of a large number of grid points (xg)  makes sense.
    if( verbose){
      print( dim( xg))
      print( dim( Z))
    }
#    
# Next call is fragile as it assumes the predictSE method includes
# a  Z argument  (even if it will often just be NULL )
      out0<- predictSE(object, xg, Z=Z, ...)
# coerce back into image format      
      out <-  as.surface( xg, out0)
    #
    # if extrapolate is FALSE set all values outside convex hull to NA
    if (!extrap) {
        if( is.null( object$x)){
          stop("need and x matrix in object")
        }
        if (is.na(chull.mask)) {
            chull.mask <- unique.matrix(object$x[, xy])
        }
        out$z[!in.poly(xg[, xy], xp = chull.mask, convex.hull = TRUE)] <- NA
    }
    #
    return(out)
}




