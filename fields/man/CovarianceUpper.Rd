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
\name{CovarianceUpper}
\alias{ExponentialUpper}
\title{
Evaluate covariance over upper triangle of distance matrix
}
\description{
Evaluates the covariance over the upper triangle of a distance matrix 
  rather than over the entire matrix to reduce computation time.  Note 
  that the \code{chol} function only requires the upper triangle of 
  the covariance matrix to perform the Cholesky decomposition.
}
\usage{
ExponentialUpper(distMat, range = 1, alpha = 1/range)
}
\arguments{
  \item{distMat}{
The distance matrix to evaluate the covariance over.
}
  \item{range}{
Range parameter default is one. Note that the scale can also
  be specified through the "aRange" scaling argument used in
  fields covariance functions)
}
  \item{alpha}{
1/range
}
}
\value{
The covariance matrix, where only the upper triangle is calculated.
}
\author{
John Paige
}
\seealso{
\code{\link[fields]{Exponential}}
}
\examples{
set.seed(123)

#make distance matrix using the random locations
coords = matrix(runif(10), ncol=2)
distMat = rdist(coords)

#compute covariance matrix, but only over the upper triangle
upperCov = ExponentialUpper(distMat, range=.1)

print(distMat)
print(upperCov)
}
\keyword{ covariance }
