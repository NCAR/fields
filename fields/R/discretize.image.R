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
"discretize.image" <- function(x, m = 64, n = 64, 
    grid = NULL, expand = c(1+ 1e-8, 1+1e-8),
    boundary.grid = FALSE, na.rm=TRUE) {
#     
# NOTE m and n are ignored if grid is passed with $x and $y.  
    if (length(expand) == 1){ 
        expand <- rep(expand, 2)
    }
    if (is.null(grid)) {
      # if grid is not given  create a  midpoint grid first
        xr <- range(x[, 1], na.rm = na.rm)
        deltemp <- (xr[2] - xr[1]) * (expand[1] - 1) * 0.5
        gridX <- seq(xr[1] - deltemp, xr[2] + deltemp, , m)
        yr <- range(x[, 2], na.rm = na.rm)
        deltemp <- (yr[2] - yr[1]) * (expand[2] - 1) * 0.5
        gridY <- seq(yr[1] - deltemp, yr[2] + deltemp, , n)
        grid <- list(x= gridX, y=gridY)
        # convert to boundary grid if needed
        if( boundary.grid){
          grid$x<- fields.convert.grid(grid$x)
          grid$y<- fields.convert.grid(grid$y)
        }
    }
    if (!boundary.grid) {
# find cut points for boundaries assuming grid has midpoints      
        xcut <- fields.convert.grid(grid$x)
        ycut <- fields.convert.grid(grid$y)
    }
    else {
        # cut points given boundaries
        xcut <- grid$x
        ycut <- grid$y
    }
    # at this point the xcut and ycut are a boundary grid 
    # even if passed as midpoints. 
    # bin ids for each location
#    cat("call to discretize.image", fill=TRUE)
    withinGrid<- (x[,1] >= min(xcut)) & (x[,1] <= max(xcut)) &
                  (x[,2] >= min(ycut)) & (x[,2] <= max(ycut)) 
    
    ##################FIX THIS
    ###################
    index1 <- findInterval(x[, 1], xcut , 
                                left.open=FALSE, rightmost.closed=TRUE) 
    index2 <- findInterval(x[, 2], ycut ,
                                left.open=FALSE, rightmost.closed=TRUE )
    if( any(!withinGrid)){
      warning("Some locations are outside the range of the grid")
      index1[!withinGrid]<- NA
      index2[!withinGrid]<- NA
    }
    tempHist<- table( index1, index2)
    ix<- as.numeric(dimnames( tempHist)[[1]])
    iy<- as.numeric(dimnames( tempHist)[[2]])
# 2 d histogram of locations
    mBin    <- length( xcut) - 1
    nBin    <- length( ycut) - 1
    hist<- matrix( 0, mBin,nBin)
    hist[ix,iy] <- tempHist 
#  save discretized locations
    xMidpoints<- (xcut[1:mBin] + xcut[2:(mBin+1)])/2
    yMidpoints<- (ycut[1:nBin] + ycut[2:(nBin+1)])/2
    loc <- cbind( xMidpoints[ index1[withinGrid] ],
                  yMidpoints[ index2[withinGrid] ] )  
    return( list( m=mBin,n=nBin, grid=grid, index=data.frame(index1, index2),
                   withinGrid=withinGrid,
                  ix= ix, iy=iy, hist=hist, loc=loc) )
}
