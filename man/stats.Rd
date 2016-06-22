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

\name{stats}
\alias{stats}
\title{
  Calculate summary statistics 
}
\description{
Various summary statistics are calculated for different 
types of data. 
}
\usage{
stats(x, by)
}
\arguments{
\item{x}{
The data structure to compute the statistics. This can either be a 
vector, matrix (data sets are the columns), or a list (data sets are 
the components). 
}
\item{by}{
If x is a vector, an optional vector (either character or numerical) 
specifying the categories to divide x into separate data sets. 
}
}
\value{
A matrix where rows index the summary statistics and the columns 
index the separate data sets. 
}
\details{
Stats breaks x up into separate data sets and then calls describe 
to calculate the statistics. 
Statistics are found by columns for matrices, by components for a list and
by the relevent groups when a  numeric vector and a by  vector are given.
The default set of statistics are the number of 
(nonmissing) observations, mean, standard deviation, minimum, lower quartile, 
median, upper quartile, maximum, and number of missing observations.  If any 
data set is nonnumeric, missing values are returned for the statistics. 
The by argument is a useful way to calculate statistics on parts of a 
data set according to different cases. 
}
\seealso{
stats.bin, stats.bplot, describe
}
\examples{
#Statistics for 8 normal random samples: 
zork<- matrix( rnorm(200), ncol=8) 
stats(zork) 

zork<- rnorm( 200)
id<- sample( 1:8, 200, replace=TRUE)
stats( zork, by=id)
}
\keyword{univar}
% docclass is function
% Converted by Sd2Rd version 1.21.
