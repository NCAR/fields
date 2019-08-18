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

\name{fields-stuff}
\alias{fields.diagonalize2}
\alias{fields.diagonalize}
\alias{fields.duplicated.matrix}
\alias{fields.mkpoly}
\alias{fields.derivative.poly}
\alias{fields.evlpoly}
\alias{fields.evlpoly2}

\title{Fields supporting functions}
\description{
Some supporting functions that are internal to fields top level
methods. Variants of these might be found in the R base but these
have been written for cleaner code or efficiency.
}
\usage{
fields.diagonalize2(A,B, verbose=FALSE)
fields.diagonalize(A,B)
fields.duplicated.matrix(mat, digits = 8) 

fields.mkpoly(x, m = 2)

fields.derivative.poly(x, m,dcoef)

fields.evlpoly( x, coef)

fields.evlpoly2( x, coef, ptab)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{A}{
A positive definite matrix}
\item{B}{
A positive definite matrix}
\item{mat}{
Arbitrary matrix for examining rows}
\item{digits}{Number of significant digits to use for comparing
elements to determine duplciate values. }

\item{x}{
Arbitrary matrix where rows are components of a multidimensional
vector}
\item{m}{ The null space degree  -- results in a polynomial of degree
(m-1) }
\item{dcoef}{ Coefficients of a multidimensional polynomial}

\item{coef}{Polynomial coefficients.}
\item{ptab}{Table of powers of different polnomial terms.}   
\item{verbose}{If TRUE prints condition number of A+B}
}

\details{
Given two matrices A (positive definite) and B (nonnegative definite)
 \code{fields.diagonalize} and \code{fields.diagonalize2} finds the 
matrix transformation G that will convert A to a identity matrix and B 
to a diagonal matrix:

 G\^T A G= I       G\^T B G= D. 


\code{fields.diagonalize2} is not as easy to follow as \code{fields.diagonalize} but may be more stable
 and is the version used in 
the Krig engine. 

\code{fields.duplicated} finds duplicate rows in a matrix.  The digits
arguments is the number of digits that are considered in the
comparison. 
The returned value is an array of integers from 1:M where M is the
number of unique rows and duplicate rows are referenced in the same
order that they appear as the rows of \code{mat}. 

\code{fields.mkpoly}  computes the complete matrix of all monomial
terms up to degree (m-1). Each row of \code{x} is are the componets of
a vector.  (The fields function mkpoly returns the number of these
terms.) In 2 dimensions with m=3 there   6 polynomial terms  up to
quadratic ( 3-1 =2) order and will be returned as the matrix:

cbind(  1 , x[,1], x[,2], x[,1]**2, x[,1]*x[,2], x[,2]**2  )

This function is used for the fixed effects polynomial  or spatial
drift used in spatial estimating functions Krig, Tps and mKrig. 
The matrix ptab is a table of the powers in each term for each
variable and is included as an attribute to the matrix returned by
this function. 
See the \code{attr} function for extracting an attribute from an
object. 

\code{ptab} for the example above is 
\preformatted{
    [,1] [,2]
[1,]    0    0
[2,]    1    0
[3,]    0    1
[4,]    2    0
[5,]    1    1
[6,]    0    2
}

This information is used in finding derivatives of the polynomial. 

\code{fields.deriviative.poly} finds the partial derivative matrix of
a multidimensional polynomial of degree (m-1) at different vector
values and with coefficients \code{dcoef}. 
This function has been orgainzed to be a clean utility for the
predicting the derivative of the estimated function from Krig or
mKrig. 
Within the fields context
the polynomial itself would be evaluated as 
     fields.mkpoly( x,m)\%*\%dcoef.
If x has  d columns ( also the dimension of the polynomial) and n rows
the partial derivatives of this polynomial at the locations x can be
organized in a nXd matrix.  This is the object returned by ths
function.   

\code{evlpoly} and \code{evlpoly2} are FORTRAN based functions for 
evaluating univariate polynomials and multivariate polynomials.  The
table of powers (ptab) needed for evlpoly2 is the same format as that
returned my the fields.mkpoly function. 


}
\author{Doug Nychka}
\seealso{Krig, Tps, as.image, predict.Krig, predict.mKrig, 
Krig.engine.default, Wendland}
\keyword{spatial}
% at least one, from doc/KEYWORDS
