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

\name{Krig.null.function}
\alias{Krig.null.function}
\title{Default function to create fixed matrix part of spatial process model.}
\description{
Constructs a matrix of terms representing a low order polynomial 
and binds additional columns due to covariates ( the Z matrix)
}
\usage{
Krig.null.function(x, Z = NULL, drop.Z = FALSE, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Spatial locations }
  \item{Z}{ Other covariates to be associated with each location.}
  \item{drop.Z}{If TRUE only the low order polynomial part is created.  }
  \item{m}{ The polynomial order is (m-1). }
}
\value{
A matrix where the first columns are the polynomial terms and the
following columns are from Z. 

}
\details{
This function can be modified to produce a different fixed part of the 
spatial model. The arguments x, Z and drop.Z are required but other arguments
can  be passed as part of a list in null.args  in the call to Krig. 
}
\author{Doug Nychka }
\seealso{Krig}
\keyword{ spatial}% at least one, from doc/KEYWORDS
