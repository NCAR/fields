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
\name{supportsArg}
\alias{supportsArg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tests if function supports a given argument
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Tests if the given function supports the given argument.  Commonly 
  used in fields code for determining if a covariance function 
  supports precomputation of the distance matrix and evaluation of 
  the covariance matrix over only the upper triangle.
}
\usage{
supportsArg(fun=stationary.cov, arg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fun}{
The function tested for support for whether it supports the argument 
  \code{arg} as input
}
  \item{arg}{
The argument to check if \code{fun} supports using as input
}
}
\details{
Currently only \code{stationary.cov}  and \code{Exp.cov} support 
  evaluation of the covariance matrix over the upper triangle 
  (and diagonal) only via the onlyUpper argument and distance 
  matrix precomputation via the distMat argument.
}

\value{
A logical indicating whether the given function supports use of the 
  given argument
}
\author{
%%  ~~who you are~~
John Paige
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{stationary.cov}}, \code{\link{Exp.cov}}
These covariance functions have the \code{onlyUpper} option allowing 
  the user to evaluate the covariance matrix over the upper triangle and 
  diagonal only and to pass a precomputed distance matrix
}
\examples{
################
#Test covariance function to see if it supports evaluation of 
#covariance matrix over upper triangle only
################

supportsArg(Rad.cov, "distMat")
supportsArg(Rad.cov, "onlyUpper")
supportsArg(stationary.cov, "distMat")
supportsArg(stationary.cov, "onlyUpper")
supportsArg(Exp.cov, "distMat")
supportsArg(Exp.cov, "onlyUpper")
}

