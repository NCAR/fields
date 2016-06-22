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

\name{fields testing scripts}
\alias{fields.tests}
\alias{test.for.zero}

\title{
Testing fields functions
}
\description{

Some of the basic methods in fields can be tested by directly
implementing the linear algebra using matrix expressions and
other functions can be cross checked within fields. These
comparisons are done in the the R source code test files in the tests
subdirectory of fields.  The function \code{test.for.zero} is useful for
comparing the tests in a meaninful and documented way. 
}

\usage{
test.for.zero( xtest, xtrue,  tol= 1.0e-8, relative=TRUE, tag=NULL)
}

\arguments{

\item{xtest}{Vector of target values}

\item{xtrue}{Vector of reference values}

\item{tol}{Tolerance to judge whether the test passes.}

\item{relative}{If true a relative error comparison is used. (See
 details below.)}

\item{tag}{ A text string to be printed out with the test results as a
   reference}

}

\details{
IMPORTANT: If the R object \code{test.for.zero.flag} exists with any
value ( e.g. \code{test.for.zero.flag <- 1} ) then when the test fails
this function will also generate an error in addition to printing a
message. This option is added to insure
that any test scripts will generate an error when any individual test
fails.  

An example:
\preformatted{
> test.for.zero( 1:10, 1:10 + 1e-10, tag="First test")
Testing:  First test
PASSED test at tolerance  1e-08
}

\preformatted{
> test.for.zero( 1:10, 1:10 + 1e-10, tag="First test", tol=1e-12)
Testing:  First test
FAILED test value =  1.818182e-10  at tolerance  1e-12
}
\preformatted{
> test.for.zero.flag <- 1
Testing:  First test
FAILED test value =  1.818182e-10  at tolerance  1e-12
Error in test.for.zero(1:10, 1:10 + 1e-10, tag = "First test", tol = 1e-12) :
}

The scripts in the \code{tests} subdirectory are 

\describe{
\item{Krig.test.R:}{Tests basic parts of the Krig and Tps functions including replicated
and weighted observations. }

\item{Krig.se.test.R:}{Tests computations of standard errors for the
Kriging estimate.}

\item{Krig.se.grid.test.R}{Tests approximate standard errors for the
Krig function found by Monte Carlo conditional simulation.} 

\item{Krig.test.W.R}{Tests predictions and A matrix  when an off diagonal 
observation weight matrix is used.}

\item{Krig.se.W.R}{Tests standard errors when an off diagonal 
observation weight matrix is used.}

\item{spam.test.R}{Tests sparse matrix formats and linear algebra.}

\item{Wend.test.R}{Tests form for Wendland covariance family and its
use of sparse matrix formats.}

\item{diag.multiply.test.R}{Tests special (efficient) version of matrix 
multiply for diagonal matrices.}

\item{ evlpoly.test.R}{Tests evaluation of univariate and multivariate
polynomial evaluation.}

\item{mKrig.test.R}{Tests the micro Krig function with and without sparse 
matrix methods. }


}

To run the tests just attach the fields library and source the testing
file.  In the  fields source code these are in a subdirectory
"tests". Compare the output to the "XXX.Rout.save" text file. 

 \code{test.for.zero} is used to print out the result for each
individual comparison.
Failed tests are potentially  bad and are reported with a
string beginning 

"FAILED test value = ... "

If the object test.for.zero.flag exists then an error is also generated
when the test fails.


FORM OF COMPARISON:
The actual test done is the sum of absolute differnces:

test value = \code{ 
sum( abs(c(xtest) - c( xtrue) ) ) /denom}

Where \code{denom} is either  \code{ mean( abs(c(xtrue)))} for relative error
or 1.0 otherwise. 

Note the use of "c" here to stack any structure in xtest and xtrue into
a vector. 

}
 
\keyword{misc}
