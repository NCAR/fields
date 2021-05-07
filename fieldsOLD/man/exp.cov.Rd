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

\name{Covariance functions} 
\alias{Exp.cov}
\alias{Exp.simple.cov}
\alias{Rad.cov}
\alias{Rad.simple.cov}
\alias{stationary.cov}
\alias{stationary.taper.cov}
\alias{wendland.cov}
\alias{cubic.cov}
\title{
  Exponential  family, radial basis 
functions,cubic spline,  compactly  supported Wendland family and 
stationary covariances. }
\description{
Given two sets of locations these functions compute the cross covariance matrix for
some covariance families. In addition these functions can take advantage 
of spareness, implement more efficient multiplcation of the 
cross covariance by a vector or matrix and also return a marginal 
variance to be consistent with calls by the Krig function.  

\code{stationary.cov} and \code{Exp.cov} have additional arguments for 
precomputed distance matrices and for calculating only the upper triangle 
and diagonal of the output covariance matrix to save time.  Also, they 
support using the \code{rdist} function with \code{compact=TRUE} or input 
distance matrices in compact form, where only the upper triangle of the 
distance matrix is used to save time.

Note: These functions have been been renamed from the previous fields functions
using 'Exp' in place of 'exp' to avoid conflict with the generic exponential 
function (\code{exp(...)})in R.
}
\usage{
Exp.cov(x1, x2=NULL, theta = 1, p=1, distMat = NA, 
  C = NA, marginal = FALSE, onlyUpper=FALSE)

Exp.simple.cov(x1, x2, theta =1, C=NA,marginal=FALSE)

Rad.cov(x1, x2, p = 1, m=NA, with.log = TRUE, with.constant = TRUE, 
               C=NA,marginal=FALSE, derivative=0)

cubic.cov(x1, x2, theta = 1, C=NA, marginal=FALSE) 

Rad.simple.cov(x1, x2, p=1, with.log = TRUE, with.constant = TRUE, 
               C = NA, marginal=FALSE)

stationary.cov(x1, x2=NULL, Covariance = "Exponential", Distance = "rdist", 
  Dist.args = NULL, theta = 1, V = NULL, C = NA, marginal = FALSE, 
  derivative = 0, distMat = NA, onlyUpper = FALSE, ...)

stationary.taper.cov(x1, x2, Covariance="Exponential", 
           Taper="Wendland", 
           Dist.args=NULL, Taper.args=NULL, 
           theta=1.0,V=NULL, C=NA, marginal=FALSE,
           spam.format=TRUE,verbose=FALSE,...)

wendland.cov(x1, x2, theta = 1, V=NULL, k = 2, C = NA, 
             marginal =FALSE,Dist.args = list(method = "euclidean"), 
             spam.format = TRUE, derivative = 0, verbose=FALSE)
}
\arguments{
\item{x1}{ Matrix of first set of locations where each row gives the
coordinates of a particular point.}

\item{x2}{ Matrix of second set of locations where each row gives the
coordinatesof a particular point. If this is missing x1 is used. }

\item{theta}{ Range (or scale) parameter. This should be a scalar (use
the V argument for other scaling options). Any distance calculated for
a covariance function is divided by theta before the covariance function
is evaluated.}
 
\item{V}{ A matrix that describes the inverse linear transformation of
the coordinates before distances are found.  In R code this
transformation is: \code{x1 \%*\% t(solve(V))} Default is NULL, that
is the transformation is just dividing distance by the scalar value
\code{theta}. See Details below.  If one has a vector of "theta's"
that are the scaling for each coordinate then just express this as
\code{V = diag(theta)} in the call to this function.}

\item{C}{ A vector with the same length as the number of rows of x2.
If specified the covariance matrix will be multiplied by this vector.}

\item{marginal}{If TRUE returns just the diagonal elements of the
covariance matrix using the \code{x1} locations. In this case this is
just 1.0. The marginal argument will trivial for this function is a
required argument and capability for all covariance functions used
with Krig.}

\item{p}{ Exponent in the exponential covariance family.  p=1 gives an
exponential and p=2 gives a Gaussian.  Default is the exponential
form.  For the radial basis function this is the exponent applied to
the distance between locations.}
\item{m}{For the radial basis function p = 2m-d, with d being the dimension of the
locations, is the exponent applied to
the distance between locations. (m is a common way of parametrizing this exponent.)}

\item{with.constant}{ If TRUE includes complicated constant for radial
 basis functions.  See the function \code{radbad.constant} for more
 details. The default is TRUE, include the constant. Without the usual
 constant the lambda used here will differ by a constant from spline
 estimators ( e.g.  cubic smoothing splines) that use the
 constant. Also a negative value for the constant may be necessary to
 make the radial basis positive definite as opposed to negative
 definite. }

\item{with.log}{If TRUE include a log term for even dimensions.  This
is needed to be a thin plate spline of integer order. }

\item{Covariance}{Character string that is the name of the covariance
shape function for the distance between locations. Choices in fields
are \code{Exponential}, \code{Matern}}

\item{Distance}{Character string that is the name of the distance
function to use. Choices in fields are \code{rdist},
\code{rdist.earth}}

\item{Taper}{Character string that is the name of the taper function
to use. Choices in fields are listed in help(taper).}

\item{Dist.args}{ A list of optional arguments to pass to the Distance
function.}

\item{Taper.args}{ A list of optional arguments to pass to the Taper
function. \code{theta} should always be the name for the range (or
scale) paremeter.}

\item{spam.format}{If TRUE returns matrix in sparse matrix format
implemented in the spam package. If FALSE just returns a full
matrix. }

\item{k}{The order of the Wendland covariance function.  See help on
 Wendland.}

 \item{derivative}{ If nonzero evaluates the partials of the
covariance function at locations x1. This must be used with  the "C" option
and is mainly called from within a predict function.  The partial
derivative is taken with respect to \code{x1}.  }

\item{verbose}{If TRUE prints out some useful information for
debugging.}

\item{distMat}{
If the distance matrix between \code{x1} and \code{x2} has already been 
computed, it can be passed via this argument so it won't need to be 
recomputed.
}
\item{onlyUpper}{
For internal use only, not meant to be set by the user.  Automatically 
set to \code{TRUE} by \code{mKrigMLEJoint} or \code{mKrigMLEGrid} if 
\code{lambda.profile} is set to \code{TRUE}, but set to \code{FALSE} 
for the final parameter fit so output is compatible with rest of 
\code{fields}.

If \code{TRUE}, only the upper triangle and diagonal of the covariance 
matrix is computed to save time (although if a non-compact distance 
matrix is used, the onlyUpper argument is set to FALSE).  If \code{FALSE}, 
the entire covariance matrix is computed.  In general, it should 
only be set to \code{TRUE} for \code{mKrigMLEJoint} and \code{mKrigMLEGrid}, 
and the default is set to \code{FALSE} so it is compatible with all of
\code{fields}.
}

\item{\dots}{ Any other arguments that will be passed to the
covariance function. e.g. \code{smoothness} for the Matern.}  }

\value{ If the argument C is NULL the cross covariance matrix is
returned.  In general if nrow(x1)=m and nrow(x2)=n then the returned
matrix will be mXn.  Moreover, if x1 is equal to x2 then this is the
covariance matrix for this set of locations.


If C is a vector of length n, then returned value is the
multiplication of the cross covariance matrix with this vector.
 
} \details{ For purposes of illustration, the function
\code{Exp.cov.simple} is provided in fields as a simple example and
implements the R code discussed below.  List this function out as a
way to see the standard set of arguments that fields uses to define a
covariance function.  It can also serve as a template for creating new
covariance functions for the \code{Krig} and \code{mKrig}
functions. Also see the higher level function \code{stationary.cov} to
mix and match different covariance shapes and distance functions.

A common scaling for stationary covariances: If \code{x1} and
 \code{x2} are matrices where \code{nrow(x1)=m} and \code{nrow(x2)=n}
 then this function will return a mXn matrix where the (i,j) element
 is the covariance between the locations \code{x1[i,]} and
 \code{x2[j,]}. The exponential covariance function is computed as
 exp( -(D.ij)) where D.ij is a distance between \code{x1[i,]} and
 \code{x2[j,]} but having first been scaled by theta. Specifically if
 \code{theta} is a matrix to represent a linear transformation of the
 coordinates, then let \code{u= x1\%*\% t(solve( theta))} and \code{v=
 x2\%*\% t(solve(theta))}.  Form the mXn distance matrix with
 elements:

\code{D[i,j] = sqrt( sum( ( u[i,] - v[j,])**2 ) )}.

and the cross covariance matrix is found by \code{exp(-D)}.  The
tapered form (ignoring scaling parameters) is a matrix with i,j entry
\code{exp(-D[i,j])*T(D[i,j])}. With T being a positive definite
tapering function that is also assumed to be zero beyond 1.

Note that if theta is a scalar then this defines an isotropic
covariance function and the functional form is essentially
\code{exp(-D/theta)}.

Implementation: The function \code{r.dist} is a useful FIELDS function
that finds the cross Euclidean distance matrix (D defined above) for
two sets of locations. Thus in compact R code we have

  exp(-rdist(u, v))

Note that this function must also support two other kinds of calls:

If marginal is TRUE then just the diagonal elements are returned (in R
code \code{diag( exp(-rdist(u,u)) )}).

If C is passed then the returned value is \code{ exp(-rdist(u, v))
\%*\% C}.


Some details on particular covariance functions:

\describe{ \item{Radial basis functions (\code{Rad.cov}:}{The
functional form is Constant* rdist(u, v)**p for odd dimensions and
Constant* rdist(u,v)**p * log( rdist(u,v) ) For an m th order thin plate
spline in d dimensions p= 2*m-d and must be positive. The constant,
depending on m and d, is coded in the fields function
\code{radbas.constant}. This form is only a generalized covariance
function -- it is only positive definite when restricted to linear
subspace.  See \code{Rad.simple.cov} for a coding of the radial basis
functions in R code.}

\item{Stationary covariance \code{stationary.cov}:}{Here the
computation is to apply the function Covariance to the distances found
by the Distance function.  For example

\code{Exp.cov(x1,x2, theta=MyTheta)}

and

\code{stationary.cov( x1,x2, theta=MyTheta, Distance= "rdist",
Covariance="Exponential")}

are the same. This also the same as

\code{stationary.cov( x1,x2, theta=MyTheta, Distance= "rdist",
Covariance="Matern",smoothness=.5)}.  }

\item{Stationary tapered covariance \code{stationary.taper.cov}:}{The
resulting cross covariance is the direct or Shure product of the
tapering function and the covariance. In R code given location
matrices, \code{x1} and \code{x2} and using Euclidean distance.
 
\code{Covariance(rdist( x1, x2)/theta)*Taper( rdist( x1,
x2)/Taper.args$theta)}

By convention, the \code{Taper} function is assumed to be identically
zero outside the interval [0,1]. Some efficiency is introduced within
the function to search for pairs of locations that are nonzero with
respect to the Taper. This is done by the SPAM function
\code{nearest.dist}.  This search may find more nonzero pairs than
dimensioned internally and SPAM will try to increase the space. One
can also reset the SPAM options to avoid these warnings.  For
spam.format TRUE the multiplication with the \code{C} argument is done
with the spam sparse multiplication routines through the "overloading"
of the \code{\%*\%} operator.  }

}

About the FORTRAN: The actual function \code{Exp.cov} and 
\code{Rad.cov} call FORTRAN to 
make the evaluation more efficient this is especially important when the 
C argument is supplied. So unfortunately the actual production code in 
Exp.cov is not as crisp as the R code sketched above. See  
\code{Rad.simple.cov} for a R coding of the radial basis functions.  

}

\seealso{
 Krig, rdist, rdist.earth, gauss.cov, Exp.image.cov, Exponential, Matern, 
Wendland.cov, mKrig} 

\examples{
# exponential covariance matrix ( marginal variance =1) for the ozone
#locations 
out<- Exp.cov( ChicagoO3$x, theta=100)

# out is a 20X20 matrix

out2<- Exp.cov( ChicagoO3$x[6:20,],ChicagoO3$x[1:2,], theta=100)
# out2 is 15X2 matrix 

# Kriging fit where the nugget variance is found by GCV 
# Matern covariance shape with range of 100.
# 

fit<- Krig( ChicagoO3$x, ChicagoO3$y, Covariance="Matern", theta=100,smoothness=2)

data( ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]
# Omit the NAs
good<- !is.na( y)
x<- x[good,]
y<- y[good]


# example of calling the taper version directly 
# Note that default covariance is exponential and default taper is 
# Wendland (k=2).

stationary.taper.cov( x[1:3,],x[1:10,] , theta=1.5, Taper.args= list(k=2,theta=2.0,
                       dimension=2) )-> temp
# temp is now a tapered 3X10 cross covariance matrix in sparse format. 

 is.spam( temp)  # evaluates to TRUE

# should be identical to
# the direct matrix product

 temp2<- Exp.cov( x[1:3,],x[1:10,], theta=1.5) * Wendland(rdist(x[1:3,],x[1:10,]), 
                      theta= 2.0, k=2, dimension=2)
 test.for.zero(  as.matrix(temp), temp2)

# Testing that the "V" option works as advertized ...
x1<- x[1:20,]
x2<- x[1:10,]

V<- matrix( c(2,1,0,4), 2,2)
Vi<- solve( V)

u1<- t(Vi\%*\% t(x1))
u2<- t(Vi\%*\% t(x2))

look<- exp(-1*rdist(u1,u2))
look2<- stationary.cov( x1,x2, V= V)
test.for.zero( look, look2)


# Here is an example of how the cross covariance multiply works
# and lots of options on the arguments


 Ctest<- rnorm(10)
 
 temp<- stationary.cov( x,x[1:10,], C= Ctest, 
        Covariance= "Wendland", 
            k=2, dimension=2, theta=1.5 )

# do multiply explicitly

 temp2<- stationary.cov( x,x[1:10,],
        Covariance= "Wendland",
            k=2, dimension=2, theta=1.5 )\%*\% Ctest

 test.for.zero( temp, temp2)


# use the tapered stationary version 
# cov.args is part of the argument list passed to stationary.taper.cov
# within Krig. 
# This example needs the spam package.
# 

\dontrun{

Krig(x,y, cov.function = "stationary.taper.cov", theta=1.5,
      cov.args= list(Taper.args= list(k=2, dimension=2,theta=2.0) )
           ) -> out2 
# NOTE: Wendland is the default taper here. 
}

# BTW  this is very similar to 
\dontrun{
 Krig(x,y, theta= 1.5)-> out
}



}
\keyword{spatial}
% docclass is function

