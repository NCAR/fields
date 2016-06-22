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

\name{US}
\alias{US}
\title{
  Plot of the US with state boundaries 
}
\description{
Plots quickly, medium resolution outlines of the US with the states and 
bodies of water. A simple wrapper for the map function from the maps package. 
}
\usage{
US( ...)
}
\arguments{
\item{\dots}{
These are the arguments that are passed to the map function
from the maps package.
}
}
\details{
The older version of this function (fields < 6.7.2) used the FIELDS dataset US.dat for the
  coordinates. Currenty this has been switched to use the maps package.    
}
\seealso{
world
}
\examples{
# Draw map in device color # 3 
US( col=3) 
}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
