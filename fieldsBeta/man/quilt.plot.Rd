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

\name{quilt.plot}
\alias{quilt.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Useful plots for visualizing irregular spatial data. }
\description{
Given a vector of z values associated with 2-d locations this 
function produces an 
image-like plot where the locations are discretized to a grid  and the z 
values are coded as a color level from a color scale. }
\usage{
quilt.plot(x, y, z, nx = 64, ny = 64,  grid = NULL, 
           add.legend=TRUE,add=FALSE, nlevel=64, 
           col = tim.colors(nlevel),
           nrow=NULL, ncol=NULL,FUN =
                 NULL, plot=TRUE, na.rm=FALSE, ...)
		 
bubblePlot= function( x, y, z, col=viridis(256),
                      horizontal = FALSE,
                      legend.shrink = 0.9, 
                      legend.width = 1.2, 
                      legend.mar = ifelse(horizontal, 3.1, 5.1),
                      size=1.0, add=FALSE,
                      ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of the x coordinates of the locations -or- a
     a 2 column matrix of the x-y coordinates. }
  \item{y}{A vector of the y coordinates -or- if the locations are passed 
in x the z vector  }
  \item{z}{Values of the variable to be plotted.}
  \item{nlevel}{Number of color levels.}
  \item{nx}{Number of grid boxes in x if a grid is not specified.}
  \item{ny}{Number of grid boxes in y. }
  \item{nrow}{Depreciated, same as nx.}
  \item{ncol}{Depreciated same as ny. }
  \item{grid}{A grid in the form of a \code{grid list}.  }
  \item{add.legend}{If TRUE a legend color strip is added}
  \item{add}{If FALSE add to existing plot.}
  \item{col}{Color scale for the image, the default is tim.colors --
    a pleasing spectrum for quilt.plot and the  modern and versitle
     viridis for bubblePlot }
    \item{plot}{If FALSE just returns the image object instead of plotting it.}
    \item{FUN}{The function to apply to values that are common to a grid box. The default is to find the mean. (see \code{as.image}).}
    \item{na.rm}{If FALSE NAs are not removed from zand so a grid box
    even one of these values may be an NA. (See details below.)}
  \item{\dots}{ for quilt.plot  arguments to be passed to the image.plot function and for bubblePlot arguments to be passed to the plot function.
  }
}
\details{

\strong{quilt.plot} 
This function combines the discretization to an image by the function 
\code{as.image} and is then graphed by \code{image.plot}.
By default, locations that fall into the same grid box will have their
z values averaged. This also means that observations that are NA will
result in the grid box average also being NA and can produce unexpected
results because the NA patterns can dominate the figure. If you are
unsure of the effect  try \code{na.rm = TRUE} for a comparison.

A similar function exists in the lattice package and produces spiffy
looking plots. The advantage of this fields version is that it uses
the
standard R graphics functions and is written in R code. Also, the
aggregation to average values for z values in the same grid box allows
for different choices of grids. If two locations are very close,
separating them could result in very small boxes. 

Legend placement is never completely automatic. Place the 
legend independently for more control, perhaps using \code{image.plot}
in tandem with \code{split.screen} or enlarging the plot margin 
See \code{help(image.plot)} for examples of this function and these 
strategies.

\strong{bubblePlot} Why was this function was written for fields?
Certainly ggplot has many options for this kind of figure. 
To quote Tim Hoar a gifted data scientist and software engineer at NCAR:
"because we could".   It is a crisp
implementation of this type of plot using lower level fields functions.
The user may choose  simply to use the source code as the basis for
a more
detailed function. However, this is also a quick plot to introduce in teaching.

This function is experimental in the sense it explores setting
out a plotting region in advance of the actual plotting using R base graphics.
See the functions  \code{\link{setupLegend}} and  \code{\link{addLegend}}
for more details. Other graphics approaches in R such as ggplot determine
the plotting regions and layout based on having the entire figure
specification at hand. Although this a comprehensive solution it also seems
overkill to just add a lone color bar to annotate a plot. Moreover, and the graphics lower level functions
to add the color bar legend are already available from the image.plot function.
}

\author{D.Nychka}
\seealso{ as.image, image.plot, lattice, persp, drape.plot,  } 
\examples{

data( ozone2)
# plot 16 day of ozone data set

quilt.plot( ozone2$lon.lat, ozone2$y[16,])
US( add=TRUE, col="grey", lwd=2)

bubblePlot( ozone2$lon.lat, ozone2$y[16,])
US( add=TRUE, col="grey", lwd=2)



### adding a common legend strip "by hand"
## and a custom color table

coltab<- two.colors( 256, middle="grey50" )

par( oma=c( 0,0,0,5)) # save some room for the legend
set.panel(2,2)
zr<- range( ozone2$y, na.rm=TRUE)

for( k in 1:4){
quilt.plot( ozone2$lon.lat, ozone2$y[15+k,], add.legend=FALSE,
 zlim=zr, col=coltab, nx=40, ny=40)
US( add=TRUE)
}
par( oma=c(0,0,0,1))
image.plot(zlim=zr,legend.only=TRUE, col=coltab)
# may have to adjust number of spaces in oma to make this work.
  }
\keyword{hplot}% at least one, from doc/KEYWORDS
