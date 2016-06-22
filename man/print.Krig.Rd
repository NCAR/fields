%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2016
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

\name{print.Krig}
\alias{print.Krig}
\title{
  Print kriging fit results. 
}
\description{
Prints the results from a fitting a spatial process estimate (Krig)
}
\usage{
\method{print}{Krig}(x,digits=4,...)
}
\arguments{
\item{x}{
Object from Krig function. 
}
\item{digits}{
Number of significant digits in printed output. Default is 4. 
}
\item{\dots}{ Other arguments to print.}
}
\value{
Selected summary results from Krig. 
}
\seealso{
print, summary.Krig, Krig 
}
\examples{
fit<- Krig(ChicagoO3$x,ChicagoO3$y, theta=100) 
print(fit) # print the summary 
fit # this will work too 
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
