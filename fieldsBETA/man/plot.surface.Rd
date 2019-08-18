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

\name{plot.surface}
\alias{plot.surface}
\title{
  Plots a surface 
}
\description{
Plots a surface object in several different ways to give 3-d
information e.g. a contour plots, perspective plots.  
}
\usage{
\method{plot}{surface}(x, main = NULL, type = "C", zlab = NULL, xlab = NULL,
    ylab = NULL, levels = NULL, zlim = NULL, graphics.reset = NULL,
     labcex = 0.6, add.legend=TRUE, ...)
}
\arguments{
\item{x}{
A surface object. At the minimum a list with components x,y and z 
in the same form as the input list for the standard contour, persp
or image functions. This can also be an object from predictSurface. 
}
\item{main}{
Title for plot. 

}
\item{type}{
type="p" for a perspective/drape plot (see drape.plot), 
type="I" for an image plot with a legend
strip (see image.plot), type="c" draws a contour plot,  type="C" is the "I" option but with contours lines
added. type="b" gives both  "p" and "C" as a 2X1 panel }
\item{zlab}{
z-axes label 
}
\item{xlab}{
x-axes label 
}
\item{ylab}{
y-axes labels 
}
\item{levels}{
Vector of levels to be passed to contour function. 
}
\item{graphics.reset}{
Reset to original graphics parameters after function plotting. 
Default is to reset if type ="b" but not for the single plot options.  
}
\item{zlim}{
Sets z limits on perspective plot. }
\item{labcex}{
Label sizes for axis labeling etc. 
}
\item{add.legend}{ If TRUE adds a legend to the draped perspective plot}
\item{\dots}{
Other graphical parameters that are passed along to either drape.persp 
or image.plot 
}

}
\seealso{
surface, predictSurface, as.surface, drape.plot, image.plot  
}
\examples{
x<- seq( -2,2,,80)
y<- seq( -2,2,,80)
# a lazy way to create some test image
z<- outer( x,y, "+")

# create basic image/surface object
obj<- list(x=x, y=y,z=z)

# basic contour plot
# note how graphical parameters appropriate to contour are passed
plot.surface( obj, type="c", col="red")

# using a fields function to fit a surface and evaluate as surface object.
fit<- Tps( BD[,1:4], BD$lnya) # fit surface to data 
# surface of variables 2 and 3 holding 1 and 4 fixed at their median levels
 out.p<-predictSurface(fit, xy=c(2,3))  

 plot.surface(out.p) # surface plot  

}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
