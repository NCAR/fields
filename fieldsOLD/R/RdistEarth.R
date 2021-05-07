# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2019
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

RdistEarth <- function(x1, x2=NULL, miles=TRUE, R=NULL){
    stopifnot(is.numeric(x1), is.matrix(x1), ncol(x1)==2,
              is.null(x2) || (is.numeric(x2) && is.matrix(x2) && ncol(x2)==2),
              (!is.null(R) && is.numeric(R)) || is.logical(miles))
    if(is.null(R)) 
        R <- if(miles[1]) 3963.34 else 6378.388
    if(is.null(x2)){
        ans <- numeric(nrow(x1)^2)
        .Call("distMatHaversin", p1=x1, radius=R, ans=ans)
        attr(ans, "dim") <- c(nrow(x1), nrow(x1))
        return(ans)
    }
    ans <- numeric(nrow(x1)*nrow(x2))
    .Call("distMatHaversin2", p1=x1, p1=x2, radius=R, ans=ans)
    attr(ans, "dim") <- c(nrow(x1), nrow(x2))
    ans
}

