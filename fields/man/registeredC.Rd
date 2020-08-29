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
\name{registeringCode}
\alias{addToDiagC}
\alias{ExponentialUpperC}
\alias{compactToMatC}
\alias{multebC}
\alias{multwendlandg}
\alias{mltdrb}
\alias{RdistC}
\alias{distMatHaversin}
\alias{distMatHaversin2}
%
\docType{data}
\title{Information objects that register
C and FORTRAN functions.
}
\description{
These are objects of class \code{CallRoutine} or \code{FortranRoutine} and also \code{\link{NativeSymbolInfo}}
They provide information for compiledfunctions called with \code{.Call},
or \code{.Fortran}.
Ordinarily one would not need to consult these and they are used to make
the search among dynamically loaded libraries ( in particular the fields library) have less ambiguity and
also be faster. These are created when the package/library is loaded 
are have their definitions from the compliation of  \code{init.c} in the
package source (src) directory. 
}
%ExponentialUpperC.Rd	compactToMatC.Rd	multebC.Rd
%RdistC.Rd		mltdrb.Rd		multwendlandg.Rd

%%\usage{
%data(addToDiagC)

%data(ExponentialUpperC)

%data(compactToMatC)

%data(multebC)

%data(multwendlandg)

%data(mltdrb)

%data(RdistC)
%}
\format{
  The format is a list with components:
  \describe{
  \item{name}{The (registration ?) name of the C function.}
  \item{address}{See \link{NativeSymbolInfo}. }
  \item{dll}{Dynamically linked library information.}
  \item{numParameters}{Number of calling arguments in function.}
  }
  }
\details{
 \describe{
\item{addToDiagC}{ adds diagonal elements to a matrix. See code{mKrig}.}

\item{ExponentialUpperC}{Fills in upper triangle of a matrix with the exponential covariance function. See \code{ExponentialUpper}}

\item{compactToMatC}{ Converts compact format to full matrix format. See \code{compactToMat}.}

\item{multebC}{Mulitplies a vector/matrix with an exponential covariance function. See \code{exp.cov}} 

\item{multwendlandg}{This has been mysteriously included but it is not a function! }

\item{mltdrb}{Evaluates the derivatives of thin plate sline radial basis functions. See \code{rad.cov}. }

\item{RdistC}{ Euclidean distance function between sets of coordinates.
See \code{rdist}.}

\item{distMatHaversin}{ Used in \code{RdistEarth}.}

\item{distMatHaversin2}{ Used in \code{RdistEarth}.}

}

See \code{package_native_routine_registration_skeleton} for the 
utility used to create these data objects.

It is not clear why these routines have been flagged as needing
documentation while other routines have not. 

}
\references{
For background on registering C, C++ and Fortran functions see 5.4 of
Writing R Extensions. 
See \url{https://r.789695.n4.nabble.com/Registration-of-native-routines-td4728874.html} for additional dicsussion of code registration. 
}
\examples{
print(addToDiagC)

}
\keyword{datasets}
