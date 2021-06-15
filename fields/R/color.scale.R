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
color.scale <- function(z, col = tim.colors, NC=256,
    zlim = NULL, transparent.color = "white", eps = 1e-08) {
#
# converts real values to a color scale of NC values.
# role of eps is to prevent values exactly at the end of the range from being
# missed
#    
# convert factor or character z to integers
    colFactor<- FALSE
    levelsZ<- NULL
    
    if( is.character(z)){
        z<- as.factor( z)
    }
    
    if( is.factor( z)){
        levelsZ<- levels(z)
        z<- as.integer(z) 
        coFactor<-TRUE
        NC<- length( levelsZ)
    }

    if (is.null(zlim)) {
        zlim <- range(z, na.rm = TRUE)
    }
    z[(z < zlim[1]) | (z > zlim[2])] <- NA
    
    if( is.function(col)){
        col<- col(NC)
    }
    
    NC <- length(col)
    span<- zlim[2]-zlim[1]
    # expand breaks slightly to include obs on the boundaries. 
    breaks <- seq(zlim[1] - span*eps,
                  zlim[2] + span*eps , 
                    length.out= NC + 1)
    # the magic of R ...
    icolor <- cut(c(z), breaks)@.Data
    # returned values is a vector of character hex strings encoding the colors.
     colorMap<- ifelse(is.na(icolor),
                         transparent.color, col[icolor])
# add zlim and col table as attributes so one can add
#  a legend and know what one is doing.     
    attr( colorMap,"zlim")<- zlim
    attr(colorMap,"col")<- col
    attr(colorMap,"levelsZ")<- levelsZ
# potential hooks for colorBar functions
    colorMapInfo<- list( col=col, zlim=zlim)
# .fieldsColorMapInfo <<- colorMapInfo
    return( colorMap)
}

