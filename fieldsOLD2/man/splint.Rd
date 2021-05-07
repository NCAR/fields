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

\name{splint}
\alias{splint}
\title{
  Cubic spline interpolation 
}
\description{
A fast, FORTRAN based function for cubic spline interpolation. 
}
\usage{
splint(x, y, xgrid, wt = NULL, derivative = 0, lam = 0, df =
                 NA, lambda = NULL, nx = NULL, digits = 8)
}
\arguments{
\item{x}{
The x values that define the curve or a two column matrix of  
x and y values. 
}
\item{y}{
The y values that are paired with the x's. 
}
\item{xgrid}{
The grid to evaluate the fitted cubic interpolating curve. 
}
\item{derivative}{
Indicates whether the function or a a first or second derivative 
should be evaluated. 
}
\item{wt}{Weights for different obsrevations in the scale of reciprocal 
variance.}
\item{lam}{ Value for smoothing parameter. Default value is zero giving 
interpolation.}
\item{lambda}{Same as \code{lam} just to make this easier to remember.}
\item{df}{ Effective degrees of freedom. Default is to use lambda =0 or a  
df equal to the number of observations.}
\item{nx}{If not NULL this should be the number of points
to evaluate on an equally spaced grid in the
range of \code{x}}
\item{digits}{Number of significant digits uused to determine what is a replicate x value.}

}
\value{
A vector consisting of the spline evaluated at the grid values in \code{xgrid}. 
}
\details{
Fits a piecewise interpolating or smoothing cubic 
polynomial to the x and y values. 
This code is designed to be fast but does not many options in  
\code{sreg} or other more statistical implementations.  
To make the solution well posed the
the second and third derivatives are set to zero at the limits of the  x 
values. Extrapolation outside the range of the x 
values will be a linear function. 

It is assumed that there are no repeated x values; use sreg followed by
predict if you do have replicated data. 

}
\section{References}{
See Additive Models by Hastie and Tibshriani. 
}
\seealso{
sreg,  Tps  
}
\examples{
x<- seq( 0, 120,,200)

# an interpolation
splint(rat.diet$t, rat.diet$trt,x )-> y

plot( rat.diet$t, rat.diet$trt)
lines( x,y)
#( this is weird and not appropriate!)

# the following two smooths should be the same

splint( rat.diet$t, rat.diet$con,x, df= 7)-> y1

# sreg function has more flexibility than splint but will
# be slower for larger data sets. 

sreg( rat.diet$t, rat.diet$con, df= 7)-> obj
predict(obj, x)-> y2 

# in fact predict.sreg interpolates the predicted values using splint!

# the two predicted lines (should) coincide
lines( x,y1, col="red",lwd=2)
lines(x,y2, col="blue", lty=2,lwd=2)
 
}
\keyword{smooth}
% docclass is function
% Converted by Sd2Rd version 1.21.
