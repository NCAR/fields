%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2018
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

\name{The Engines:}
\alias{Krig.engine.default}
\alias{Krig.engine.knots}
\alias{Krig.engine.fixed}
\alias{Krig.coef}
\alias{Krig.check.xY}
\alias{Krig.cor.Y}
\alias{Krig.transform.xY}
\alias{Krig.make.u}
\alias{Krig.make.W}
\alias{Krig.make.Wi}
\alias{\%d*\%}
\title{ Basic linear algebra utilities and other 
    computations supporting the Krig function. }
\description{
These are internal functions to Krig that compute the basic matrix
decompositions or solve the linear systems needed to evaluate the 
Krig/Tps estimate. Others listed below do some simple housekeeping and 
formatting. Typically they are called from within Krig but can also be
used directly if passed a Krig object list. 
}
\usage{
Krig.engine.default(out, verbose = FALSE)
Krig.engine.knots(out, verbose = FALSE)
Krig.engine.fixed( out, verbose=FALSE, lambda=NA)

Krig.coef(out, lambda = out$lambda, y = NULL, yM = NULL, verbose = FALSE) 
Krig.make.u(out, y = NULL, yM = NULL, verbose = FALSE) 
Krig.check.xY(x, Y,Z, weights, na.rm, verbose = FALSE) 
Krig.cor.Y(obj, verbose = FALSE) 
Krig.transform.xY(obj, knots, verbose = FALSE)

Krig.make.W( out, verbose=FALSE)
Krig.make.Wi ( out, verbose=FALSE)

}
\arguments{

 \item{out}{ A complete or partial Krig object. If partial it must have
all the information accumulated to this calling point within the Krig
function.  }

\item{obj}{Same as \code{out}. }

\item{verbose}{If TRUE prints out intermediate results for
debugging.}

\item{lambda}{Value of smoothing parameter "hard wired" into decompositions. 
Default is NA, i.e.  use the value in \code{out\$lambda}. }


\item{y}{New y vector for recomputing coefficients. OR for \%d*\% a
vector or matrix. }

\item{yM}{New y vector for recomputing coefficients but 
the values have already been collapsed into replicate group means.}

\item{Y}{raw data Y vector}

\item{x}{raw x matrix of spatial locations  
OR 
In the case of \%d*\%, y is either a matrix or a vector. As a vector, y, is
interpreted to be the elements of a digaonal matrix. 
}

\item{weights}{ Raw \code{weights} vector passed to Krig}
\item{Z}{ Raw vector or matrix of additional covariates.}

\item{na.rm}{ NA action logical values passed to Krig}

\item{knots}{Raw \code{knots} matrix  passed to Krig} 

}

\details{

ENGINES:

The engines are the 
code modules that handle the basic linear algebra needed to 
computed the estimated curve or surface coefficients. 
All the engine work on the data that has been reduced to unique 
locations and possibly replicate group means with the weights 
adjusted accordingly. All information needed for the decomposition are 
components in the Krig object passed to these functions. 

 \code{Krig.engine.default} finds the decompositions for a Universal
Kriging estimator. by simultaneously diagonalizing the linear system
system for the coefficients of the estimator. The main advantage of this
form is that it is fairly stable numerically, even with ill-conditioned
covariance matrices with lambda > 0. (i.e.  provided there is a "nugget"
or measure measurement error. Also the eigendecomposition allows for
rapid evaluation of the likelihood, GCV and coefficients for new data
vectors under different values of the smoothing parameter, lambda. 


 \code{Krig.engine.knots} finds the decompositions in the case that the
covariance is evaluated at arbitrary locations possibly different than
the data locations (called knots). The intent of these decompositions is
to facilitate the evaluation at different values for lambda.  There will
be computational savings when the number of knots is less than the
number of unique locations. (But the knots are as densely distributed as
the structure in the underlying spatial process.) This function call
fields.diagonalize, a function that computes the matrix and eigenvalues
that simultaneous diagonalize a nonnegative definite and a positive
definite matrix. These decompositions also facilitate multiple
evaluations of the likelihood and GCV functions in estimating a
smoothing parameter and also multiple solutions for different y vectors. 

 \code{Krig.engine.fixed} are specific decomposition based on the Cholesky 
factorization assuming that the smoothing parameter is fixed. This 
is the only case that works in the sparse matrix.
Both knots and the full set of locations can be handled by this case. 
The difference between the "knots" engine above is that only a single value
of lambda is considered in the fixed engine. 


OTHER FUNCTIONS:

\code{Krig.coef} Computes the "c" and "d" coefficients to represent the 
estimated curve. These coefficients are used by the predict functions for 
evaluations. Krig.coef can be used outside of the call to Krig to 
recompute the fit with different Y values and possibly with different
lambda values. If new y values are not passed to this function then the yM
vector in the Krig object is used. The internal function 
\code{Krig.ynew} sorts out the logic of what to do and use based on the 
passed arguments. 

\code{Krig.make.u} Computes the "u" vector, a transformation of the collapsed
observations that allows for rapid evaluation of the GCV function and 
prediction. This only makes sense when the decomposition is WBW or DR, i.e. 
an eigen decomposition. If the decompostion is the Cholesky based then this
function returns NA for the u component in the list.  

\code{Krig.check.xY} Checks for removes missing values (NAs).

\code{Krig.cor.Y} Standardizes the data vector Y based on a correlation model. 

 \code{Krig.transform.xY} Finds all replicates and collapse to unique
locations and mean response and pooled variances and weights. These are
the xM, yM and weightsM used in the engines. Also scales the x locations
and the knots according to the transformation. 


 \code{Krig.make.W} and \code{Krig.make.Wi} These functions create an
off-diagonal weight matrix and its symmetric square root or the inverse
of the weight matrix based on the information passed to Krig.  If
\code{out$nondiag} is TRUE W is constructed based on a call to the passed
function wght.function along with additional arguments.  If this flag is
FALSE then W is just \code{diag(out$weightsM)} and the square root and inverse
are computed directly. 


\code{\%d*\%} Is a simple way to implement efficient diagonal
multiplications.  x\%d*\%y is interpreted to mean diag(x)\%*\% y
if x is a vector. If x is a matrix then this becomes the same as the usual 
matrix multiplication. 
}

\section{Returned Values}{

ENGINES:

The returned value is a list with the matrix decompositions and 
other information. These are incorporated into the complete Krig object. 

Common to all engines:
\describe{
   \item{decomp}{Type of decomposition}
   \item{nt}{dimension of T matrix}
   \item{np}{number of knots}
}

\code{Krig.engine.default}:

\describe{
\item{u}{Transformed data using eigenvectors.}
\item{D}{Eigenvalues}
\item{G}{Reduced and weighted matrix of the eigenvectors}
\item{qr.T}{QR decomposition of fixed regression matrix}
\item{V}{The eigenvectors}
}

\code{Krig.engine.knots}:

\describe{
  \item{u}{A transformed vector that is based on the data vector.}
  \item{D}{Eigenvalues of decomposition}
  \item{G}{Matrix from diagonalization}
  \item{qr.T}{QR decomposition of the matrix for the fixed component. 
 i.e. sqrt( Wm)\%*\%T}
   \item{pure.ss}{pure error sums of squares including both the
    variance from replicates and also the sums of squared residuals
    from fitting the full knot model with lambda=0 to the replicate means. }
}

\code{Krig.engine.fixed}:

\describe{
\item{d}{estimated coefficients for the fixed part of model}
\item{c}{estimated coefficients for the basis functions derived from the 
                 covariance function.}
}

Using all data locations

\describe{
\item{qr.VT}{QR decomposition of the inverse Cholesky factor times the 
T matrix. }
\item{MC}{Cholesky factor}
}

Using knot locations
\describe{
\item{qr.Treg}{QR decomposition of regression matrix modified by the 
estimate of the nonparametric ( or spatial) component.} 
\item{lambda.fixed}{Value of lambda used in the decompositions}
}

OTHER FUNCTIONS:

\code{Krig.coef}
\describe{
\item{yM}{Y values as replicate group means}
\item{shat.rep}{Sample standard deviation of replicates}
\item{shat.pure.error}{Same as shat.rep}
\item{pure.ss}{Pure error sums of squares based on replicates}
\item{c}{The "c" basis coefficients associated with the covariance
or radial basis functions.}
\item{d}{The "d" regression type coefficients that are from the fixed part of the model
or the linear null space.}
\item{u}{When the default decomposition is used the data vector transformed by the orthogonal matrices. This facilitates evaluating the GCV function 
at different values of the smoothing parameter.}
}              
\code{Krig.make.W}
\describe{
\item{W}{The weight matrix}
\item{W2}{ Symmetric square root of weight matrix}
}

\code{Krig.make.Wi}
\describe{
\item{ Wi}{The inverse weight matrix}
\item{W2i}{ Symmetric square root of inverse weight matrix}
}
 

}

\author{Doug Nychka }

\seealso{  \code{\link{Krig}}, \code{\link{Tps}} }
\examples{

Krig( ChicagoO3$x, ChicagoO3$y, theta=100)-> out

Krig.engine.default( out)-> stuff

# compare "stuff" to components in out$matrices

look1<- Krig.coef( out)
look1$c
# compare to out$c

look2<- Krig.coef( out, yM = ChicagoO3$y)
look2$c
# better be the same even though we pass as new data!

}
\keyword{ spatial }
