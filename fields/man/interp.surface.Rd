%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2018
%# University Corporation for Atmospheric Research (UCAR)
%# Contact: Douglas Nychka, nychka@mines.edu,
%# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
%#
%# This program is free software; you can redistribute it and/or modify
%# it under the terms of the GNU General Public License as published by
%# the Free Software Foundation; either version 2 of the License, or
%# (at your option) any later version.
%# This program is distributed in the hope that it will be useful,
%# but WITHOUT ANY WARRANTY; without even the implied warranty of
%# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%# GNU General Public License for more details.
%#
%# You should have received a copy of the GNU General Public License
%# along with the R software environment if not, write to the Free Software
%# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%# or see http://www.r-project.org/Licenses/GPL-2    

\name{interp.surface}
\alias{interp.surface}
\alias{interp.surface.grid}
\title{
  Fast bilinear interpolator from a grid.
}
\description{
Uses bilinear weights to interpolate values on a rectangular
grid to arbitrary locations or to another grid.
}
\usage{
interp.surface(obj, loc)
interp.surface.grid(obj, grid.list)

}
\arguments{
\item{obj}{
A list with components x,y, and z in the same style as used by contour,
persp, image etc. x and y are the X and Y grid values and z is a matrix
with the corresponding values of the surface
}
\item{loc}{
A matrix of (irregular) locations to interpolate. First column of loc
isthe X coordinates and second is the Y's. 
}
\item{grid.list}{ A list with components x and y 
describing the grid to interpolate. The grids do not need to be equally spaced.}

}
\value{
 An vector of interpolated values.  NA are returned for
regions of the obj\$z that are NA and also for locations outside of the
range of the parent grid. 
}
\details{
Here is a brief explanation of
the interpolation:  Suppose that the location, (locx, locy)  lies in
between the first two grid points in both x an y. That is  locx is between
x1 and x2 and
locy is between y1 and y2.  Let ex= (l1-x1)/(x2-x1) ey= (l2-y1)/(y2-y1).
The
interpolant is

( 1-ex)(1-ey)*z11 + (1- ex)(ey)*z12 + ( ex)(1-ey)*z21 + ( ex)(ey)*z22  

Where the z's are the corresponding elements of the Z matrix. 

Note that bilinear interpolation can produce some artifacts related to
the grid and not reproduce higher behavior in the surface. For, example
the extrema of the interpolated surface will always be at the parent
grid locations.  There is nothing special about about interpolating to
another grid, this function just includes a \code{for} loop over one
dimension and a call to the function for irregular locations. It was
included in fields for convenience. since the grid format is so common. 

See also the akima package for fast interpolation from irrgeular locations. 
Many thanks to  Jean-Olivier Irisson for making this code more efficient and 
concise. 

}
\seealso{
image.smooth, as.surface, as.image, image.plot, krig.image,Tps
}
\examples{
#
# evaluate an image at a finer grid
# 

data( lennon)
# create an example in the right list format like image or contour
obj<- list( x= 1:20, y=1:20, z= lennon[ 201:220, 201:220])

set.seed( 123)
# lots of random points
N<- 500
loc<- cbind( runif(N)*20, runif(N)*20)
z.new<- interp.surface( obj, loc)
# compare the image with bilinear interpolation at scattered points
set.panel(2,2)
image.plot( obj)
quilt.plot( loc, z.new) 


# sample at 100X100 equally spaced points on a grid

grid.list<- list( x= seq( 1,20,,100), y=  seq( 1,20,,100))

interp.surface.grid( obj, grid.list)-> look

# take a look
set.panel(2,2)
image.plot( obj)
image.plot( look)

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
