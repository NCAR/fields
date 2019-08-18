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

\name{tim.colors}
\alias{tim.colors}
\alias{larry.colors}
\alias{two.colors}
\alias{designer.colors}

\alias{color.scale}
\alias{snow.colors}
\alias{fieldsPlotColors}
\title{ Some useful color tables for images and tools to handle them. }
\description{
Several color scales useful for image plots: a pleasing rainbow style
color table patterned after that used in Matlab by Tim Hoar and also
some simple color interpolation schemes between two or more
colors. There is also a function that converts between colors and a
real valued vector.
}
\usage{
tim.colors(n = 64, alpha=1.0)

larry.colors()

snow.colors(n=256, alpha=1)

two.colors(n=256, start="darkgreen", end="red", middle="white",
alpha=1.0)

designer.colors( n=256, col= c("darkgreen", "white", "darkred"), x=
                              seq(0,1,, length(col)) ,alpha=1.0)

color.scale( z, col=tim.colors(256), zlim =NULL,
transparent.color="white", eps= 1e-8)

fieldsPlotColors( col,...)

}
\arguments{
\item{alpha}{The transparency of the color -- 1.0 is opaque and 0 is
             transparent.  This is useful for overlays of color and
             still being able to view the graphics that is covered. }

\item{n}{ Number of color levels. The setting \code{n}=64 is the
orignal definition.}

\item{start}{Starting color for lowest values in color scale}

\item{end}{ Ending color.}

\item{middle}{Color scale passes through this color at halfway}

\item{col}{A list of colors (names or hex values) to interpolate}

\item{x}{Positions of colors on a [0,1] scale. Default is to assume
that the x values are equally spacesd from 0 to 1.}

\item{z}{Real vector to encode in a color table.}

\item{zlim}{Range to use for color scale. Default is the
\code{range(z)} inflated by 1- eps and 1+eps.}

\item{transparent.color}{Color value to use for NA's or values outside
\code{zlim}}

\item{eps}{A small inflation of the range to avoid boundary values of
\code{z} being coded as NAs}
\item{\dots}{Additional plotting arguments to code{image.plot}}

}

\details{
The color in R can be represented as three vectors in RGB coordinates
and these coordinates are interpolated separately using a cubic spline
to give color values that intermediate to the specified colors.

 Ask Tim Hoar about \code{tim.colors}! He is a Mattlab black belt and
this is his favorite scale in that system.  \code{two.colors} is
really about three different colors. For other colors try
\code{fields.color.picker} to view possible choices.
\code{start="darkgreen", end="azure4"} are the options used to get a
nice color scale for rendering aerial photos of ski trails.  (See
\url{https://github.com/dnychka/MJProject}.)  \code{larry.colors}
is a 13 color palette used by Larry McDaniel (retired software engineer from NCAR) and is particularly
useful for visualizing fields of climate variables.

\code{snow.colors} is the scale used by Will Klieber's team for
visualizing snow cover from remotely sensed data products. See the commented code for the script as to how how this was formed from an orignal raw 256 level scale. Note the that first color in this table is grey and is desigend to represent the minimum value of the range ( e.g. 0). If the image in in percent snow cover then \code{zlim=c(0,100)}
would make sense as a range to fit grey pixels to zero and white to 100 percent.

 \code{designer.color} is the master function for the otther scales. 
 It can be useful if one wants to customize the color
table to match quantiles of a distribution.  e.g. if the median of the
data is at .3 with respect to the range then set \code{x} equal to
c(0,.3,1) and specify three colors to provide a transtion that matches
the median value. In fields language this function interpolates
between a set of colors at locations x. While you can be creative
about these colors just using another color scale as the basis is
easy. For example

\code{designer.color( 256, rainbow(4), x= c( 0,.2,.8,1.0))}

leaves the choice of the colors to Dr. R after a thunderstorm.
See also colorBrewer to choose sequences of colors that form a good palette. 

\code{color.scale} assigns colors to a numerical vector in the same way as 
the \code{image} function. This is useful to kept the assigment of colors consistent
across several vectors by specifiying a common \code{zlim} range.  

\code{plotColorScale} A simple function to plot a vector of colors to examine their values. 
}

\value{

A vector giving the colors in a hexadecimal format, two extra hex
digits are added for the alpha channel.

}
\seealso{ topo.colors, terrain.colors, image.plot, quilt.plot, grey.scale,
fields.color.picker }
\examples{

tim.colors(10) 
# returns an array of 10 character strings encoding colors in hex format

# e.g. (red, green,  blue) values of   (16,255, 239)
#   translates to "#10FFEF" 
# rgb( 16/255, 255/255, 239/255, alpha=.5)
#   gives   "#10FFEF80"  note extra "alpha channel"

# view some color table choices
set.panel( 4,1)
fieldsPlotColors( tim.colors())
title("tim.colors")
fieldsPlotColors( larry.colors())
title("larry.colors")
fieldsPlotColors( two.colors())
title("two.colors")
fieldsPlotColors( snow.colors())
title("snow.colors")

# a bubble plot with some transparency for overlapping dots
set.seed(123)
loc<- matrix( rnorm( 200), 100,2)
Z<- loc[,1] + loc[,2]
colorMap<- color.scale( Z, col=tim.colors(10, alpha=.8))
par( mar=c(5,5,5,5)) # extra room on right for color bar
plot( loc, col=colorMap, pch=16, cex=2)
#  add a color scale
 image.plot(legend.only=TRUE, zlim=range( Z), col=tim.colors(10))

# using tranparency without alpha the image plot would cover points

obj<- list( x= 1:8, y=1:10, z= outer( 1:8, 1:10, "+") )
plot( 1:10,1:10)

image(obj, col=two.colors(alpha=.5), add=TRUE)

coltab<- designer.colors(col=c("blue", "grey", "green"),
                   x= c( 0,.3,1) )
		   

image( obj, col= coltab )

# peg colors at some desired quantiles  of data.
# NOTE need 0 and 1 for the color scale to make sense
x<- quantile( c(obj$z), c(0,.25,.5,.75,1.0) )
# scale these to [0,1]
zr<- range( c(obj$z))
x<- (x-zr[1])/ (zr[2] - zr[1])  

coltab<- designer.colors(256,rainbow(5), x)
image( obj$z, col= coltab )
# see image.plot for adding all kinds of legends


set.panel()

}
\keyword{ aplot}
