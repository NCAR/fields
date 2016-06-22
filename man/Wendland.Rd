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

\name{Wendland}
\alias{Wendland}
\alias{Wendland.beta}
\alias{Wendland2.2}
\alias{fields.D}
\alias{fields.pochdown}
\alias{fields.pochup}
\alias{wendland.eval}
%\alias{Wendland.father}
%\alias{Wendland.mother}
%\alias{wendland.basis}


\title{Wendland family of covariance functions and supporting 
numerical functions}

\description{
 Computes the compactly supported, stationatry Wendland covariance
function as a function ofdistance. This family is useful for creating
sparse covariance matrices. 

}

\usage{
Wendland(d, theta = 1, dimension, k,derivative=0, phi=NA)

Wendland2.2(d, theta=1)
Wendland.beta(n,k)
wendland.eval(r, n, k, derivative = 0)
fields.pochup(q, k)
fields.pochdown(q, k)
fields.D(f,name,order = 1)

%Wendland.father(x, theta = 1, dimension = 1, k=3) 
%Wendland.mother(x, theta = 1, dimension = 1, k=3) 
%wendland.basis(x1, x2, theta = 1, V=NULL,
%    k = 3, C = NA, Dist.args = list(method = "euclidean"), 
%    spam.format = TRUE, verbose = FALSE, flavor=0)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{d}{Distances between locations. Or for wendland.coef the dimension of
the locations.}
\item{theta}{Scale for distances. This is the same as the
range parameter.}
\item{dimension}{Dimension of the locations}
\item{n}{Dimension for computing Wendland polynomial coefficients}
\item{k}{Order of covariance function.}
\item{derivative}{Indicates derivative of covariance function}
\item{phi}{Depreciated argument will give stop if not an NA.
 (Formerly the scale factor to multiply the function. Equivalent to the
marginal variance or sill if viewed as a covariance function.) }
\item{r}{ Real value in [0,1] to evaluate Wendland function.}
\item{q}{Order of Pochhammer symbol}
\item{f}{Numerical expression to differentiate.}
\item{name}{Variable with which to take derivative.}
\item{order}{Order of derivative.}

%\item{x}{Argument for one dimensional basis function}
%\item{x1}{Two dimensional locations to evaluate the basis functions}
%\item{x2}{Two dimensional centers that define the basis}
%\item{C}{Vector or matrix of coefficients to multiply with basis functions}
%\item{Dist.args}{Arguments to distance function}
%\item{V}{See explantion in help for \code{Exp.cov}}
%\item{spam.format}{If TRUE return result in sparse format}
%\item{verbose}{If TRUE prints out intermediate steps for debugging.}
%\item{flavor}{Takes values 0:3. Controls type of tensor product:
%father*father=0, father*mother =1, mother*father = 2, mother*mother =3 }

}

\details{
 This is the basic function applied to distances and called by the
\code{wendland.cov} function. It can also be used as the Covariance or 
Taper specifications in the more general 
stationary.cov and station.taper.cov functions.
The proofs and construction of the Wendland family of positive definite functions can be found in the work of Wendland(1995).
( H. Wendland. Piecewise polynomial , positive definite and compactly supported radial functions of minimal degree. AICM 4(1995), pp 389-396.)

The Wendland covariance function is a positive 
polynomial on [0,theta] and zero beyond theta. It is further normalized in these fields functions to be 1 at 0. The parameter \code{k} detemines the smoothness of the covariance at zero. The additional parameter \code{n} or \code{dimension} is needed because the property of 
positive definitness for radial functions depends on the dimension being considered.  

The polynomial terms of the Wenland function.
 are computed recursively based on the values of \code{k}
and \code{dimension} in the function \code{wendland.eval}. The matrix of
coefficients found by \code{Wendland.beta} is used to weight each polynomial term and follows Wendland's original construction of these functions. The recursive definition of the Wendland coefficients depends on Pochhammer symbols akin to binomial coefficients:

\code{fields.pochup(q, k)}
calculates the Pochhammer symbol for rising factorial q(q+1)(q+2)...(q+k-1)

and 

\code{fields.pochdown(q, k)}
calculates the Pochhammer symbol for falling factorial q(q-1)(q-2)...(q-k+1). 

Derivatives are found symbolically using  a recursive modification of the base function \code{D} (\code{fields.D}) and then evaluated numerically based on the polynomial form.


A specific example of the Wendland family is \code{Wendland2.2} (k=2, dimension=2). This is included mainly for testing but the explicit formula may also be enlightening. 
}

\value{
A vector of the covariances or its derivative. 

}
\author{Doug Nychka, Ling Shen}
\seealso{ wendland.cov, stationary.taper.cov}
\examples{

dt<- seq( 0,1.5,, 200)

y<- Wendland( dt, k=2, dimension=2)

plot( dt, y, type="l")

# should agree with 

y.test<- Wendland2.2( dt)
points( dt, y.test)

# second derivative
plot( dt, Wendland( dt, k=4, dimension=2, derivative=2), type="l")

# a radial basis function using the Wendland  the  "knot" is at (.25,.25)
gl<- list( x= seq( -1,1,,60), y = seq( -1,1,,60) )

bigD<- rdist( make.surface.grid( gl), matrix( c(.25,.25), nrow=1))
RBF<- matrix(Wendland( bigD, k=2, dimension=2), 60,60)

# perspective with some useful settings for shading.
persp( gl$x, gl$y, RBF, theta=30, phi=20, shade=.3, border=NA, col="grey90")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
