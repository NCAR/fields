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

\name{ribbon.plot}
\alias{ribbon.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
 \title{ Adds to an existing plot, a ribbon of color, based on values from
a color scale, along a sequence of line segments.}
 \description{
Given a series of 2-d points and values at these segments, 
the function colors the segments according to a color scale and the
segment values. This is essentially an image plot restricted to line segments.
}
\usage{

ribbon.plot(x,y,z,zlim=NULL, col=tim.colors(256),
               transparent.color="white",...)

}

\arguments{
  \item{x}{x locations of line segments}
  \item{y}{y locations of line segments}
  \item{z}{ Values associated with each segment.}

  \item{zlim}{Range for z values to determine color scale. }

  \item{col}{Color table used for strip. Default is our favorite 
tim.colors being a scale from a dark blue to dark red.}

\item{transparent.color}{Color used for missing values. Default is that
missing values make the ribbon transparent.}
 
 \item{\dots}{Optional graphical arguments that are passed to the
   \code{segment} plotting function. A favorite is lwd to make a broad
 ribbon. }

}

\details{
Besides possible 2-d applications, 
this function is useful to annotate a curve on a surface using colors. 
The values mapped to acolor scheme could indicate a feature 
other than the height of the surface. 
For example, this function could indicate the slope of the surface. 

 }
\author{Doug Nychka}

\seealso{ image.plot, arrow.plot, add.image, colorbar.plot}

\examples{
plot( c(-1.5,1.5),c(-1.5,1.5), type="n")
temp<- list( x= seq( -1,1,,40), y= seq( -1,1,,40))
temp$z <- outer( temp$x, temp$y, "+")
contour( temp, add=TRUE)

t<- seq( 0,.5,,50)
y<- sin( 2*pi*t)
x<- cos( pi*t)
z<- x + y

ribbon.plot( x,y,z, lwd=10)

persp(  temp, phi=15, shade=.8, col="grey")-> pm
trans3d( x,y,z,pm)-> uv
ribbon.plot( uv$x, uv$y, z**2,lwd=5)

}

\keyword{ hplot }% at least one, from doc/KEYWORDS
