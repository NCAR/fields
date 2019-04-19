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

\name{flame}
\alias{flame}
\title{
  Response surface experiment ionizing a reagent 
}
\description{
The characteristics of an ionizing flame are varied with the intent of 
maximizing the intensity of emitted light for lithuim in 
solution. Areas outside of the measurements are where the mixture may 
explode! Note that the optimum is close to the boundary. Source of data is 
from a master's level lab experiment in 
analytical chemistry from Chuck Boss's course at NCSU. 
<s-section name= "DATA DESCRIPTION">
This is list with the following components
}
\arguments{
\item{x}{
x is a 2 column matrix with the different Fuel and oxygen flow rates for 
the burner. 
}
\item{y}{
y is the response. The intensity of light at a particular 
wavelength indicative of Lithium ions. 
}
}
\keyword{datasets}
% docclass is data
% Converted by Sd2Rd version 1.21.
