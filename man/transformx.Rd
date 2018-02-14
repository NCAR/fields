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

\name{transformx}
\alias{transformx}
\title{
  Linear transformation   
}
\description{
Linear transformation of each column of a matrix. There are several 
choices of the type of centering and scaling.  
}
\usage{
transformx (x, scale.type = "unit.sd", x.center, x.scale)
}
\arguments{
\item{x}{
Matrix with columns to be transformed.
}
\item{scale.type}{
Type of transformation the default is "unit.sd": subtract the mean and 
divide by the standard deviation. 
Other  choices are "unscaled" (do nothing),  "range" (transform to 
[0,1]),"user" (subtract a supplied location and divide by a scale).
}
\item{x.center}{
A vector of centering values to subtract from each column. 
}
\item{x.scale}{
A vector of scaling values to subtract from each column.
}
}
\value{
A matrix whose columns have between transformed.
This matrix also has the attributes: scale.type, x.center and y.center 
with the transformation information. 
}
\details{
After deciding what the centering and scaling values should be for each 
column of x, this function just calls the standard utility scale. 
This function was created partly to attach the transformation information 
as attributes to the transformed matrix. It is used in Krig, cover.design, 
krig.image etc. to transform the independent variables.   
}
\seealso{
scale
}
\examples{
#
newx<-transformx( ChicagoO3$x, scale.type="range")
}
\keyword{manip}
% docclass is function
% Converted by Sd2Rd version 1.21.
