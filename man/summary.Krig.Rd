%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2017
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

\name{summary.Krig}
\alias{summary.Krig}
\alias{summary.spatialProcess}
\title{
  Summary for Krig or spatialProcess estimated models. 
}
\description{
Creates a list of summary results including estimates for  the nugget
variance (sigma) and the smoothing
parameter (lambda). This list is usually printed using  a "print.summary"
function for nice formatting.
}
\usage{
\method{summary}{Krig}(object, digits=4,...)
\method{summary}{spatialProcess}(object, digits=4,...)

}
\arguments{
\item{object}{
A Krig or spatialProcess object. 
}
\item{digits}{
Number of significant digits in summary. 
}
\item{\dots}{Other arguments to summary}
}
\value{
Gives a summary of the Krig object. The components 
include the function call, number of observations, effective degrees  
of freedom, residual degrees of freedom, root mean squared error,  
R-squared and adjusted R-squared, log10(lambda), cost, GCV minimum and  
a summary of the residuals.  
}
\details{
This function is a method for the generic function summary for class
Krig. The results are formatted and printed using print.summary.Krig.
}
\seealso{
Krig, summary, print.summary.Krig  
}
\examples{
fit<- Krig(ChicagoO3$x, ChicagoO3$y, theta=100)  
summary(fit)                            # summary of fit 
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
