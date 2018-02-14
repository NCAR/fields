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

\name{grid list}
\alias{grid list}
\alias{grid.list}
\alias{fields.x.to.grid}
\alias{parse.grid.list}
\alias{fields.convert.grid}
\alias{discretize.image}
\alias{make.surface.grid}
\alias{unrollZGrid}

\title{
Some simple functions for working with gridded data and 
the grid format (grid.list) used in fields.
}
\description{
    
The object grid.list refers to a list that contains
information for evaluating a function on a 2-dimensional
grid of points. If a function has more than two
independent variables then one also needs to specify the
constant levels for the variables that are not being
varied.  This format is used in several places in fields
for functions that evaluate function estimates and plot
surfaces. These functions provide some default conversions
among information and the gird.list. The function 
\code{discretize.image} is a useful tool for "registering"
irregular 2-d points to a grid. 
}

\usage{
parse.grid.list( grid.list, order.variables="xy")
fields.x.to.grid(x,nx=80, ny=80, xy=c(1,2))
fields.convert.grid( midpoint.grid )
discretize.image(x, m = 64, n = 64, grid = NULL, 
                 expand = c(1 + 1e-08, 1 + 1e-08),
                 boundary.grid = FALSE, na.rm = TRUE)
make.surface.grid( grid.list)
unrollZGrid( grid.list, ZGrid)
}
\arguments{
\item{grid.list}{ No surprises here -- a grid list! These can be unequally 
spaced.}
\item{order.variables}{ If "xy" the x variable will be subsequently plotted
as the horizontal variable. If "yx" the x variable will be on the vertical 
axis.}
\item{x}{ A matrix of independent variables such as the locations of 
observations given to Krig.}
\item{nx}{Number of grid points for x variable.}
\item{ny}{Number of grid points for y variable.}
\item{m}{Number of grid points for x variable.}
\item{n}{Number of grid points for y variable.}
\item{na.rm}{Remove missing values if TRUE}
\item{xy}{The column positions that locate the x and y variables for the 
grid.} 
\item{grid}{ A grid list!}
\item{expand}{ 
A scalar or two column vector that will expand the
grid beyond the range of the observations.}

\item{midpoint.grid}{ Grid midpoints to convert to grid boundaries.}

\item{boundary.grid}{ If TRUE interpret grid points as boundaries  of 
grid boxes. If FALSE interpret as the midpoints of the boxes. }

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
}
\details{The form of a grid.list is  

\code{list( var.name1= what1  , var.name2=what2 , ... var.nameN=what3)}

Here var.names are the names of the independent variables. 
The what options describe what should be done with this variable when  
generating the grid. These should either an increasing sequence of points
or a single vaules. Obviously there should be only be two variables with 
sequences to define a grid for a surface.

Most of time the gridding  sequences  are equally  
spaced and are easily generated using the \code{seq} function. Also throughout fields
the grid points are typically the midpoints of the grid rather the grid box
boundaries. However, these functions can handle unequally spaced grids and the
logical boundary.grid can indicate a grid being the box boundaries. 

The variables in the list components are assumed  
to be in the same order as they appear in the data matrix. 

A useful function that expands the grid from the grid.list description into 
a full set of locations is \code{make.surface.grid} and  is 
just a wrapper around the R base function \code{expand.grid}. A typical operation is to go from a grid.list to the set of grid locations. Evaluate a
fucntion at these lcoations and then reformat this as an image for plotting. 
Here is how to do this cleanly:
\preformatted{
grid.list<- list( x= 1:10, y=1:15)
xg<- make.surface.grid(grid.list)
# look at a surface dependin on xg locations
z<- xg[,1] + 2*xg[,2]
out<- list( x=grid.list$x, y= grid.list$y, z=matrix( z, nrow=10, ncol=15))
# now for example
image.plot( out)
}

The key here is that \code{xg} and \code{matrix} both  organize the grid in the 
same order. 

Some fields internal functions that support interpreting grid list format are:

\code{fields.x.to.grid}: 
 Takes an "x" matrix of locations or independent variables and creates a
reasonable grid list. This is used to evaluate predicted surfaces when a
grid list is not explicited given to predictSurface. The variables
(i.e. columns of x) that are not part of the grid are set to the median
values.  The x grid values are \code{nx} equally spaced points in the
range \code{x[, xy[1]]}. The y grid values are \code{ny} equally spaced
points in the range \code{x[, xy[2]]}.


\code{parse.grid.list}:
 Takes a grid list and returns the information in a more expanded list
form that is easy to use. This is used, for example, by predictSurface
to figure out what to do!


\code{fields.convert.grid}:
 Takes a vector of n values assumed to be midpoints of a grid and 
returns the n+1 boundaries. See how this is used in discretize.image 
with the cut function. This function will handle unequally spaced
grid values. 

 \code{discretize.image}: Takes a vector of locations and a 2-d grid and
figures out to which boxes they belong. The output matrix ind has the
grid locations. If boundary.grid is FALSE then the grid list (grid) is
assumed to be grid midpoints. The grid boundaries are taken to be the
point half way between these midpoints. The first and last boundaries
points are determined by extrapolating so that the first and last box
has the midpoint in its center. (See the code in fields.convert.grid for
details.) If grid is NULL then midpoints are found from m and n and the
range of the x matrix. 

\code{unrollZGrid} Checks that the ZGrid object is compatible with th e grid.list and concatenates the grid arrays into vectors. This version of the covariates are used the usual predict function. 

}

\seealso{
as.surface, predictSurface, plot.surface, surface, 
expand.grid, as.image 
}

\examples{

#Given below are some examples of grid.list objects and the results  
#when they are used with make.surface.grid. Note that  
#make.surface.grid returns a matrix that retains the grid.list  
#information as an attribute. 

grid.l<- list( 1:3, 2:5) 
make.surface.grid(grid.l)

  
grid.l <- list( 1:3, 10, 1:3) 
make.surface.grid(grid.l) 

#The next  example shows how the grid.list can be used to  
#control surface plotting and evaluation of an estimated function. 
# first create a test function  

set.seed( 124)

X<- 2*cbind( runif(30), runif(30), runif(30)) -1
  
dimnames( X)<- list(NULL, c("X1","X2","X3")) 
y<- X[,1]**2 + X[,2]**2 + exp(X[,3])   

# fit an  interpolating thin plate spline  
out<- Tps( X,y) 

grid.l<- list( X1= seq( 0,1,,20), X2=.5, X3=seq(0,1,,25)) 
surface( out, grid.list=grid.l) 
#  surface plot based on a 20X25 grid in X1 an X3  
#                       over the square [0,2] and [0,2]   
#                       holding X2 equal to 1.0. 
#


# test of discretize to make sure points on boundaries are counted right
set.seed(123)
x<- matrix( runif(200), 100,2)
look<- discretize.image( x, m=2,n=2)
xc<- seq(min(x[,1]), max(x[,1]),,5)
xc<- xc[2:4]
yc<- seq(min(x[,2]), max(x[,2]),,5)
yc<- yc[2:4]
grid <-  list( x= xc, y= yc)
look2<- discretize.image( x, m=2,n=2)

table( look$index )
table( look2$index )

# indicator image of discretized locations
look<- discretize.image( RMprecip$x, m=15, n=15)
image.plot( look$grid$x, look$grid$y,look$hist )  
# actual locations
points( RMprecip$x,col="magenta", pch=".") 



}
\keyword{misc}
% docclass is function
% Converted by Sd2Rd version 1.21.
