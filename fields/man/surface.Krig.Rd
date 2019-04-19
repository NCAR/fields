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

\name{surface.Krig}
\alias{surface.Krig}
\alias{surface.mKrig}

\title{
  Plots a surface and contours 
}
\description{
Creates different plots of the fitted surface of a Krig object. This is a
quick way to look at the fitted function over reasonable default
ranges.  
}
\usage{
\method{surface}{Krig}(
object, grid.list = NULL, extrap = FALSE,
                 graphics.reset = NULL, xlab = NULL, ylab = NULL, main
                 = NULL, zlab = NULL, zlim = NULL, levels = NULL, type
                 = "C", nx = 80, ny = 80, ...)

\method{surface}{mKrig}(
object, grid.list = NULL, extrap = FALSE,
                 graphics.reset = NULL, xlab = NULL, ylab = NULL, main
                 = NULL, zlab = NULL, zlim = NULL, levels = NULL, type
                 = "C", nx = 80, ny = 80, ...)
}
\arguments{
\item{object}{
A Krig object or an mKrig object.
}
\item{grid.list}{
A list with as many components as variables describing the surface. 
All components should have a single value except the two that give the 
grid points for evaluation. If the matrix or data frame has column names,  
these must appear in the grid list. If grid.list is missing an the
surface has just two dimensions the grid is based on the ranges of the
observed data.
 
}
\item{extrap}{
Extrapolation beyond the range of the data. If false only the
convex hull of the observations is plotted. Default is false. 
}
\item{graphics.reset}{
Reset to original graphics parameters after function plotting. 
}
\item{type}{
Type of plot as a character. "p" perspective plot (persp). "c" contour
plot (contour). "b" a two panel figure with perspective and contour
plots. "I" image plot with legend strip (image.plot). "C" image plot
with contours overlaid.  Image with  contour is the default.
  
}
\item{main}{
Title of plot}
\item{xlab}{
x axis label}
\item{ylab}{
y axis label}
\item{zlab}{
z axis label if "p" or "b" type is used.}
\item{zlim}{
Z limits passed to persp}
\item{levels}{
Contour levels passed to contour. }
\item{nx}{
 Number of grid points to evaluate surface on the horizontal 
axis (the x-axis). }

\item{ny}{
 Number of grid points to evaluate surface on the vertical  
axis (the y-axis). }

\item{\dots}{
Any other plotting options. 
}
}
\details{
This function is essentially a combination of predictSurface and
plot.surface. It may not always give a great rendition but is easy to use
for checking the fitted surface. The default of extrap=F is designed to
discourage looking at the estimated surface outside the range of
the observations.  

NOTE: that any Z covariates will b edropped and only the spatial part of
the model will be evaluated. 

}
\seealso{
\code{\link{Krig}} 
predictSurface, plot.surface, image.plot
}
\examples{
fit<- Krig(ChicagoO3$x,ChicagoO3$y, theta=30)  # krig fit 

#Image plot of surface with nice, smooth  contours and shading

surface(fit, type="C", nx=128, ny=128) 
 
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
