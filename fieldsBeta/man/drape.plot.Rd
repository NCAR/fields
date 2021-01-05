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

\name{drape.plot}
\alias{drape.color}
\alias{drape.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perspective plot draped with colors in the facets.}
\description{
Function to produce the 
usual wireframe perspective plot with the facets being filled  
with different colors. By default the colors are assigned from a
color bar based on the z values. \code{drape.color}  can be used to create
a color matrix different from the z matrix used for the wireframe.}
\usage{
drape.plot(x, y, z, z2=NULL, col = tim.colors(64), zlim = range(z, na.rm=TRUE), 
 zlim2 = NULL, add.legend = TRUE, horizontal = TRUE, theta = 30, phi = 20, 
   breaks=NA, ...)

drape.color(z, col = tim.colors(64), zlim = NULL,breaks,
   transparent.color = "white", midpoint=TRUE, eps=1e-8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ grid values for x coordinate (or if x is a list the 
components x y and z are used.)}
  \item{y}{grid values for y coordinate}
  \item{z}{A matrix of z heights }
  \item{z2}{ A matrix of z values to use for coloring facets. If NULL 
then z is used for this purpose}
  \item{col}{ A color table for the z values that will be used for draping}
  \item{zlim}{the z limits for \code{z} these are used to set up
  the scale of the persp plot. This defaults to range(z, na.rm=TRUE) as
  in persp}
  \item{zlim2}{the z limits for \code{z2} these are used to set up
  the color scale. This defaults to  }
  \item{add.legend}{ If true a color strip is added as a legend.}
  \item{horizontal}{ If true color strip is put at bottom of the 
plot, if FALSE it is placed vertically on the right side.  }
  \item{theta}{ x-y rotation angle for perspective. }
  \item{phi}{ z-angle for perspective. }
  \item{transparent.color}{ Color to use when given an NA in z} 
  \item{midpoint}{ If TRUE color scale is formed for midpoints of z
       obtained by averaging 4 corners.}
  \item{breaks}{Numerical divisions for the color scale. If the default (NA)
  is N+1 equally spaced points in the range \code{zlim} where N is the number of
colors in \code{col}. This is the argument has the same effect as used in the
\code{image} and \code{image.plot} functions.}
 \item{eps}{Amount to inflate the range (1+/- eps) to inlude points on 
     break endpoints.}
  \item{\dots}{ Other arguments that will be passed to the 
          persp function. The most common is zlim the z limits for the
	  3-d  plot and also the limits to set up the color scale. The
	  default for zlim is the range of z.}
}
\value{
\code{drape.plot}
If an assignment is made the projection matrix from persp is returned. 
This information can be used to add additional 3-d features to the plot.
See the \code{persp} help file for an example how to add additional points 
and lines using the \code{trans3d} function and also the example below. 


\code{drape.color} If dim( z) = M,N this function returns a list with components:
\item{color.index}{An (M-1)X(N-1) matrix (midpoint= TRUE) or MXN matrx (midpoint=FALSE) where each element is a text string specifying the color. }
\item{breaks}{The breaks used to assign the numerical values in z to color categories.}

}
\details{ 

The legend strip may obscure part of the plot. If so, add this as
another step using image.plot.

When using \code{drape.color} just drop the results into the
\code{col} argument of \code{persp}. Given this function there are no
surprises how the higher level \code{drape.plot} works: it calls
\code{drape.color} followed by \code{persp} and finally the legend
strip is added with \code{image.plot}.

The color scales essentially default to the ranges of the z
values. However, by specifying zlim and/or zlim2 one has more control
of how the perspective plot is scaled and the limits of the color
scale used to fill the facets.  The color assignments are done by
dividing up the zlim2 interval into equally spaced bins and adding a
very small inflation to these limits.  The mean z2 values, comprising
an (M-1)X(N-1) matrix, for each facet are discretized to the bins. The
bin numbers then become the indices used for the color scale. If zlim2
is not specified it is the range of the z2 matrix is used to generate
the ranges of the color bar. Note that this may be different than the
range of the mean facets.  If z2 is not passed then z is used in its
place and in this case the zlim2 or zlim argument can used to define
the color scale.

This kind of plot is also supported through the wireframe function in the
\code{lattice} package. The advantage of the fields version is that it uses the
standard R graphics functions -- and is written in R code.  

The drape plot is also drawn by the fields \code{surface} function with 
\code{type="P"}. 

}
\author{D. Nychka }

\seealso{ image.plot, quilt.plot,  persp, plot.surface, surface, lattice, trans3d} 
\examples{

# an obvious choice:
# Dr. R's favorite New Zealand Volcano!
data( volcano)
M<- nrow( volcano)
N<- ncol( volcano)
x<- seq( 0,1,,M)
y<- seq( 0,1,,N)

pm<- drape.plot( x,y,volcano, col=terrain.colors(128)) 

# use different range for color scale and persp plot
# setting of border omits the mesh lines

 drape.plot( x,y,volcano, col=topo.colors(128),zlim=c(0,300),
                     zlim2=c( 120,200), border=NA)

# note tranparent color for facets outside the zlim2 range

#The projection has been saved in pm
# add a point marking the summit
zsummit <- max( volcano) 
ix<- row( volcano)[volcano==zsummit]
iy <- col( volcano)[volcano==zsummit]
uv <- trans3d( x[ix], y[iy],zsummit,pm)
points( uv, col="magenta", pch="+", cex=2)

# overlay volcano wireframe with gradient in x direction. 

dz<- ( 
     volcano[1:(M-1), 1:(N-1)] - volcano[2:(M), 1:(N-1)] +
     volcano[1:(M-1), 2:(N)] - volcano[2:(M), 2:(N)]  
         )/2

# convert dz to a color scale:
  zlim<- range( c( dz), na.rm=TRUE)
  zcol<-drape.color( dz, zlim =zlim, col = viridis(64) )$color.index

# with these colors 

  persp( volcano, col=zcol, theta=30, phi=20,
  border=NA,expand=.3 )

# add legend using image.plot function 
  image.plot( zlim=zlim, legend.only =TRUE, horizontal =TRUE,
            col= viridis(64))
 
}
\keyword{hplot}% at least one, from doc/KEYWORDS
