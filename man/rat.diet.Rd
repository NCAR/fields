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

\name{rat.diet}
\alias{rat.diet}
\title{
  Experiment studying an appetite supressant in rats.
}
\description{
The `rat.diet' data frame has 39 rows and 3 columns.
These are data from a study of an appetite supressant given to young rats.
The suppressant was removed from the treatment group at around 60 days.
The responses are the median food intake and  each group had approximately 
10 animals. 
}
\usage{
data(rat.diet)
}
\format{
This data frame contains the following columns:

\describe{
\item{t}{
Time in days}
\item{con}{
Median food intake of the control group }
\item{trt}{
Median food intake of the treatment group}
}

}
\keyword{datasets}
% docclass is data
% Converted by Sd2Rd version 1.21.
