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

\name{add.image}
\alias{add.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Adds an image to an existing plot.}
\description{
Adds an image to an existing plot. Simple arguments control the 
location and size. 
}
\usage{
add.image(xpos, ypos, z, adj.x = 0.5, adj.y = 0.5, 
image.width = 0.15, image.height = NULL, col = tim.colors(256), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xpos}{X position of  image in user coordinates }
  \item{ypos}{ Y position of  image in user coordinates }
  \item{z}{ Matrix of intensities comprising the image. }
 \item{adj.x}{ Location of image relative to x coordinate. 
 Most
common values are .5 (centered), 0 (right side of image at x) and  1 
(left side  of
image at x). These are the same conventions that are used for \code{adj} 
in
positioning text.}
  \item{adj.y}{Location of image relative to y coordinate.
Same rules as \code{adj.x}
}
  \item{image.width}{ Width of image as a fraction of the 
plotting region in horizontal direction. }
  \item{image.height}{ Height of image as a fraction of the
plotting region in horizontal direction. If NULL height is scaled to 
make image pixels square.}
  \item{col}{ Color table for image. Default is tim.colors.}
  \item{\dots}{Any other plotting arguments that are passed to the image 
function  }
}
\seealso{ image.plot, colorbar.plot, image, tim.colors }
\examples{
plot( 1:10, 1:10, type="n")
data( lennon)

add.image( 5,4,lennon, col=grey( (0:256)/256))
# reference lines 
xline( 5, col=2)
yline( 4,col=2) 

#
# add lennon right in the corner beyond the plotting region
# 

par(new=TRUE, plt=c(0,1,0,1), mar=c(0,0,0,0), usr=c(0,1,0,1))
add.image( 0,0, lennon, adj.x=0, adj.y=0) 

}
\keyword{ hplot }% at least one, from doc/KEYWORDS
