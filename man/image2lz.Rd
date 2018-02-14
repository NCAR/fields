%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2018
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
\name{image2lz}
\alias{image2lz}
\alias{crop.image}
\alias{in.poly}
\alias{in.poly.grid}
\alias{half.image}
\alias{get.rectangle}
\alias{average.image}
\alias{which.max.matrix}
\alias{which.max.image}
 \title{Some simple functions for subsetting images} 
\description{These function help in subsetting images or reducing its
size by averaging adjecent cells.}
\usage{

crop.image(obj, loc=NULL,...)
which.max.matrix(z)
which.max.image(obj)
get.rectangle()
average.image(obj, Q=2)
half.image(obj)
in.poly( xd, xp, convex.hull=FALSE, inflation=1e-07)
in.poly.grid( grid.list,xp, convex.hull=FALSE, inflation=1e-07)

}

\arguments{
  \item{obj}{A list in image format with the usual x,y  defining the
grid and z a matrix of image values.}

 \item{loc}{A 2 column matrix of locations within the image region
that define the subset. If not specified then the image is plotted and
the rectangle can be specified interactively.}
\item{Q}{Number of pixels to average.}
 \item{xd}{ A 2 column matrix of locations that are the points to check 
for being inside a polygon.}
 \item{xp}{ A 2 column matrix of locations that are vertices of a 
polygon. The last point is assumed to be connected to the first.}
 \item{convex.hull}{If TRUE then the convex hull of \code{xp} is used 
instead of the polygon.}
\item{grid.list}{A list with components x and y specifing the 2-d grid values.
(See help( grid.list) for more details.)}
 \item{inflation}{A small expansion factor to insure that points 
precisely on the boundaries and vertices of the convex hull are included 
as members.}
\item{z}{ A matrix of numerical values}

 \item{\dots}{ Graphics arguments passed to image.plot. This 
is only relevant when loc is NULL and the locator function is called via 
\code{get.rectangle}. }


}
\details{ If \code{loc} has more than 2 rows then the largest rectangle 
containing the locations is used. 

\describe{

\item{crop.image}{Creates a subset of the image \code{obj} by taking 
using the largest rectangle in the locations \code{loc}. This is useful 
if one needs to extract a image that is no bigger in extant than som 
edata location.  If locations are omitted the parent image is plotted 
and the locations from two mouse clicks on the image. Returned value is 
an image with appropriate \code{x,y} and \code{z} components.}

\item{get.rectangle}{Given an image plots and waits for two mouse 
clicks that are returned.}

\item{which.max.image}{Returns a list with components \code{x, y, z}
, and \code{ind} giving the
location of the maximun and value of the maximum in the image based
on the grid values and also on the indicies of the image matrix.}

\item{average.image, half.image}{Takes passed image and averages the 
pixel values and adjusts the grid to create an image that has a smaller 
number of elements. If \code{Q=2} in \code{average.image} it has the 
same effect as \code{half.image} but might be slower -- if the original 
image is mXn then half image will be an image (m/2)X(n/2). This begs the 
question what happens when m or n is odd or when (m/Q) or (n/Q) are not 
integers. In either case the largest rows or columns are dropped. (For 
large \code{Q} the function might be modified to drop about half the 
pixels at both edges.)  }

\item{in.poly, in.poly.grid}{Determines whether the points xd,yd are 
inside a polygon or outside. Return value is a logical vector with TRUE 
being inside or on boundary of polygon. The test expands the polygon 
slightly in size (on the order of single precision zero) to include 
points that are at the vertices. \code{in.poly} does not really depend 
on an image format however the grid version \code{in.poly.grid} is more 
efficient for considering the locations on a regular grid
See also \code{in.land.grid} that is hard coded to work with the
fields world map.}

}

}

\author{Doug Nychka}

\seealso{ drape.plot, image.plot,
    interp.surface, interp.surface.grid, in.land.grid}

\examples{
data(RMelevation)

# region defining Colorado Front Range

  loc<- rbind( c(-106.5, 40.8),
             c(-103.9, 37.5))

# extract elevations for just CO frontrange.
   FR<- crop.image(RMelevation, loc)
   image.plot( FR, col=terrain.colors(256))
   
   which.max.image( FR)

# average cells  4 to 1 by doing this twice!
   temp<-  half.image( RMelevation)
   temp<- half.image( temp)

# or in one step
   temp<-  average.image( RMelevation, Q=4)-> temp
   image.plot( temp, col=terrain.colors(256))

# a polygon (no special meaning entered with just locator)
x1p<- c(
 -106.2017, -104.2418, -102.9182, -102.8163, -102.8927, -103.3254, -104.7763,
 -106.5581, -108.2889, -109.1035, -109.3325, -108.7980)

x2p<- c(
  43.02978, 42.80732, 41.89727, 40.84566, 39.81427, 38.17618, 36.53810, 36.29542,
  36.90211, 38.29752, 39.45025, 41.02767)
xp<- cbind( x1p,x2p)

 image.plot( temp)
 polygon( xp[,1], xp[,2], lwd=2)

# find all grid points inside poly
 fullset<- make.surface.grid( list( x= temp$x, y= temp$y))
 ind<-  in.poly( fullset,xp)

# take a look 
 plot( fullset, pch=".")
 polygon( xp[,1], xp[,2], lwd=2)
 points( fullset[ind,], pch="o", col="red", cex=.5)

# masking out the image NA == white in the image plot
 temp$z[!ind] <- NA
 image.plot( temp)
 polygon( xp[,1], xp[,2], lwd=2)

# This is more efficient for large grids:
# because the large number of grid location ( xg above) is 
# never explicitly created.

 ind<- in.poly.grid( list( x= temp$x, y= temp$y), xp)

# now use ind in the same way as above to mask points outside of polygon

}

\keyword{ hplot }% at least one, from doc/KEYWORDS
