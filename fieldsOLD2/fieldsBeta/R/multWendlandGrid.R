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
multWendlandGrid <- function( grid.list,center, delta, coef, xy= c(1,2) ){
     xGrid<- grid.list[[xy[1]]]
     yGrid<- grid.list[[xy[2]]]
     mx<- length( xGrid)
     my<- length( yGrid)
# transform centers to correspond to integer spacing of grid:
# i.e. 1:nx and 1:ny
     dx<- (xGrid[mx] - xGrid[1]) / (mx-1)
     dy<- (yGrid[my] - yGrid[1]) / (my-1)
     centerScaled<- cbind( ((center[,1] - xGrid[1]) / dx) + 1,
                           ((center[,2] - yGrid[1]) / dy) + 1 )
     deltaX<- delta/dx
     deltaY<- delta/dy
     
     nc<- nrow( center)
     out<-.Fortran( "multWendlandG", PACKAGE="fields",
                 mx=as.integer(mx),
                 my=as.integer(my),
                 deltaX= as.double( deltaX),
                 deltaY= as.double( deltaY),                  
                 nc= as.integer(nc),
                 center=as.double(centerScaled),
                 coef=as.double(coef),
                 h= as.double(matrix(0,mx,my)),
                 flag=as.integer(1)
                )
     if( out$flag!= 0){
       stop("error in multWendlandG FORTRAN")}
     return( out$h)
   }


 