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

\name{summary.ncdf}
\alias{summary.ncdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizes a netCDF file handle}
\description{
Provides a summary of the variable names and sizes from the handle
returned from netCDF file. 
}
\usage{
\method{summary}{ncdf}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ 
The "handle" returned by the \code{read.ncdf} function from the
ncdf package.  
 }
\item{\dots}{ Other arguments to pass to this function. Currently, no 
other arguments are used. 
}
}
\details{
This function is out of place in fields but was included because often large
geophysical data sets are in netCDF format and the ncdf R package is also 
needed. To date the summary capability in the ncdf package is limited and 
this function is used to supplement it use. The function is also a a 
useful device to see how the ncdf object is structured.
}

\author{ D. Nychka }

\seealso{ ncdf}
\keyword{ IO }% at least one, from doc/KEYWORDS
