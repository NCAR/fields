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

\name{predictSurface}
\alias{predictSurface}
\alias{predictSurface.default}
\alias{predictSurface.mKrig}
\alias{predictSurface.Krig}
\alias{predictSurface.fastTps}
\alias{predictSurfaceSE}
\alias{predictSurfaceSE.default}
\alias{predict.surface}

\title{
  Evaluates a fitted function or the prediction error as a surface that is suitable for plotting with 
the image, persp, or contour functions.
}
\description{
Evaluates a a fitted model or the prediction error on a 2-D grid keeping any other variables constant.
The resulting object is suitable for use with functions for viewing 3-d
surfaces.  
}
\usage{
\method{predictSurface}{default}(object, grid.list = NULL, 
                     extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
                        xy = c(1,2),  verbose = FALSE, ...)

\method{predictSurface}{fastTps}(object, grid.list = NULL, 
       extrap = FALSE, chull.mask = NA, nx = 80, ny = 80,
       xy = c(1,2),  verbose = FALSE, ...)
       
\method{predictSurface}{Krig}(object, grid.list = NULL, extrap = FALSE, chull.mask = NA, 
nx = 80, ny = 80, xy = c(1, 2), verbose = FALSE, ZGrid = NULL, 
    drop.Z = FALSE, just.fixed=FALSE,  ...)
    
\method{predictSurface}{mKrig}(object,  ...)

\method{predictSurfaceSE}{default}(object, grid.list = NULL, extrap = FALSE, chull.mask =
                 NA, nx = 80, ny = 80, xy = c(1, 2), verbose = FALSE,
                 ZGrid = NULL, just.fixed = FALSE, ...)

\method{predict}{surface}(object,...)
}

\arguments{
\item{object}{
An object from fitting a function to data. In fields this is usually a
Krig, mKrig, or fastTps object. 
}
\item{grid.list}{
A list with as many components as variables describing the surface. 
All components should have a single value except the two that give the 
grid points for evaluation. If the matrix or data frame has column names,  
these must appear in the grid list. See the grid.list help file for more
details. If this is omitted and the fit just depends on two variables the
grid will be made from the ranges of the observed variables. 
(See the function \code{fields.x.to.grid}.)
 
}
\item{extrap}{
 Extrapolation beyond the range of the data. If \code{FALSE} (the 
default) the predictions will be restricted to the convex hull of the observed 
data or the convex hull defined from the points from the argument chull.mask. 
This function may be slightly faster if this logical is set to 
\code{TRUE} to avoid checking the grid points for membership in the 
convex hull. For more complicated masking a low level creation of a bounding 
polygon and testing for membership with \code{in.poly} may be useful. 

}

\item{chull.mask}{
Whether to restrict the fitted surface to be on a convex hull, NA's
are assigned to values outside the
convex hull. chull.mask should be a sequence of points defining a convex
hull. Default is to form the convex hull from the observations if this
argument is missing (and extrap is false).  
}

\item{nx}{
Number of grid points in X axis. }
\item{ny}{
Number of grid points in Y axis. }

\item{xy}{ A two element vector giving the positions for the "X" and "Y"
variables for the surface. The positions refer to the columns of the x 
matrix used to define the multidimensional surface. This argument is 
provided in lieu of generating the grid list. If a 4 dimensional surface
is fit to data then \code{ xy= c(2,4)} will evaluate a surface using the 
second and fourth variables with  variables 1 and 3 fixed at their median 
values. NOTE: this argument is ignored if a grid.list argument is 
passed. }
\item{drop.Z}{If TRUE the fixed part of model depending on covariates is omitted.}
\item{just.fixed}{If TRUE the nonparametric surface is omitted.}
 \item{\dots}{
Any other arguments to pass to the predict function associated with the fit
object. 
Some of the usual arguments for several of the fields fitted objects include:
\describe{

\item{ynew}{ New values of y used to reestimate the surface.}

\item{Z}{A matrix of covariates for the fixed part of model.}
}
}

\item{ZGrid}{An array  or list form of covariates to use for
	 prediction. This must match the
\code{grid.list} argument. e.g. ZGrid and grid.list describe the same
 grid. 
If ZGrid is an array then the first two indices are the x and y
 locations in the 
grid. The third index, if present, indexes the covariates. e.g. For
 evaluation on 
a 10X15 grid and with 2 covariates. \code{ dim( ZGrid) == c(10,15, 2)}.
If ZGrid is a list then the components x and y shold match those of
 grid.list and
the z component follows the shape described above for the no list 
case. 
}

 \item{verbose}{If TRUE prints out some imtermediate results for debugging.}

}

\value{
The usual list components for making contour and perspective plots
(x,y,z) along with labels for the x and y variables. For
\code{predictSurface.derivative} the component \code{z} is a three
dimensional array with \code{nx}, \code{ny}, 2.
 } 
\details{ This
function creates the right grid using the grid.list information or the
attribute in xg, calls the predict function for the object with these
points and also adding any extra arguments passed in the ... section,
and then reforms the results as a surface object (as.surface). To
determine the what parts of the prediction grid are in the convex hull
of the data the function \code{in.poly} is used. The argument
inflation in this function is used to include a small margin around
the outside of the polygon so that point on convex hull are
included. This potentially confusing modification is to prevent
excluding grid points that fall exactly on the ranges of the
data. Also note that as written there is no computational savings for
evaluting only the convex subset compared to the full grid.

\code{predictSurface.fastTps} is a specific version ( m=2, and k=2)
that can be much more efficient because it takes advantage of a low
level FORTRAN call to evaluate the Wendland covariance function. Use
\code{predictSurface} or \code{predict} for other choices of m and k.

\code{predictSurface.Krig} is designed to also include covariates for the fixed in terms of grids. Due to similarity in output and the model. \code{predictSurface.mKrig} just uses the Krig method.
 
NOTE: \code{predict.surface} has been depreciated and just prints out
a warning when called. 

 }
\seealso{
Tps, Krig, predict, grid.list, make.surface.grid, as.surface, surface, 
in.poly
}
\examples{
fit<- Tps( BD[,1:4], BD$lnya)  # fit surface to data 

# evaluate fitted surface for  first two 
# variables holding other two fixed at median values

out.p<- predictSurface(fit)
surface(out.p, type="C") 

#
# plot surface for second and fourth variables 
# on specific grid. 

glist<- list( KCL=29.77, MgCl2= seq(3,7,,25), KPO4=32.13, 
                     dNTP=seq( 250,1500,,25))

out.p<- predictSurface(fit, glist)
surface(out.p, type="C")

out.p<- predictSurfaceSE(fit, glist)
surface(out.p, type="C")
  
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
