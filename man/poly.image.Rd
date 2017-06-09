%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2017
%# University Corporation for Atmospheric Research (UCAR)
%# Contact: Douglas Nychka, nychka@ucar.edu,
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
\name{poly.image}
\alias{poly.image}
\alias{poly.image.regrid}
\title{Image plot for cells that are irregular quadrilaterals.}
\description{
 Creates an image using polygon filling based on a grid of irregular
quadrilaterals.  This function is useful for a regular grid that has
been transformed to another nonlinear or rotated coordinate system. This
situation comes up in lon-lat grids created under different map projections. 
Unlike the usual image format this function requires the grid to be
specified as two matrices x and y that given the grid x and y coordinates
explicitly for every grid point. 
}
\usage{
poly.image(x, y, z, col = tim.colors(64), breaks, transparent.color = "white", 
 midpoint = FALSE, zlim = range(z, na.rm = TRUE), 
 xlim = range(x), ylim = range(y), add = FALSE, border=NA,lwd.poly=1,...)

poly.image.regrid(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix of the x locations of the grid. }
  \item{y}{A matrix of the y locations of the grid. }
  \item{z}{Values for each grid cell. Can either be the value at the
grid points or interpreted as the midpoint of the grid cell. }
  \item{col}{ Color scale for plotting. }
  \item{breaks}{Numerical breaks to match to the colors. If missing breaks are
 equally spaced on the range \code{zlim}.}
  \item{transparent.color}{ Color to plot cells that are outside the 
range specified in the function call. }
  \item{midpoint}{ 
 Only relevant if the dimensions of x,y, and z are the same. If TRUE
the z values will be averaged and then used as the cell midpoints.  If
FALSE the x/y grid will be expanded and shifted to represent grid cells
corners. (See poly.image.regrid.) }
 
 \item{zlim}{ Plotting limits for z. }
  \item{xlim}{Plotting limits for x. }
  \item{ylim}{Plotting limits for y.}
  \item{add}{ If TRUE will add image onto current plot. }
  \item{border}{Color of the edges of the quadrilaterals, the default is
no color.}
  \item{lwd.poly}{Line width for the mesh surface. i.e. the outlines of the quadrilateral
facets. This might have to be set smaller than one if rounded corners on the facets are visible. }
  \item{\dots}{ If add is FALSE, additional graphical arguments that
 will be supplied to the plot function. }
}
\details{
 This function is straightforward except in the case when the dimensions
of x,y, and z are equal. In this case the relationship of the values to
the grid cells is ambigious and the switch midpoint gives two possible
solutions. The z values at 4 neighboring grid cells can be averaged to
estimate a new value interpreted to be at the center of the grid. This
is done when midpoint is TRUE. Alternatively the full set of z values
can be retained by redefining the grid. This is accomplisehd by finding
the midpoints of x and y grid points and adding two outside rows and
cols to complete the grid.  The new result is a new grid that is is
(M+1)X (N+1) if z is MXN. These new grid points define cells that
contain each of the original grid points as their midpoints. Of course
the advantage of this alternative is that the values of z are preserved
in the image plot; a feature that may be important for some uses. 

The function image.plot uses this function internally when image
information is passed in this format and can add a legend. In most cases
just use image.plot.

The function \code{poly.image.regrid} does a simple averaging and
extrapolation  of the grid locations to shift from midpoints to
corners.  In the interior grid corners are found by the average of the
4 closest midpoints. For the edges the corners are just extrapolated
based on the separation of nieghboring grid cells. 
}
\author{Doug Nychka }
\seealso{image.plot}
\examples{
data(RCMexample)
set.panel( 1,2)
par(pty="s")
# plot with grid modified
poly.image( RCMexample$x, RCMexample$y, RCMexample$z[,,1])

# use midpoints of z
poly.image( RCMexample$x, RCMexample$y, RCMexample$z[,,1],midpoint=TRUE)

  set.panel()
# an example with quantile breaks

 brk<- quantile(  RCMexample$z[,,1], c( 0, .9,.95,.99,1.0) )
 poly.image( RCMexample$x, RCMexample$y, RCMexample$z[,,1], breaks=brk, col=
    rainbow(4))
  

# images are very similar. 
  set.panel()
# Regridding of x and y
  l1<- poly.image.regrid( RCMexample$x)
  l2<- poly.image.regrid( RCMexample$y)

# test that this works
  i<- 1:10
  plot( l1[i,i], l2[i,i])
  points( RCMexample$x[i,i], RCMexample$y[i,i],col="red")

 

  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
