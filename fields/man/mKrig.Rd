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
\name{mKrig}
\alias{mKrig}
\alias{predict.mKrig}
\alias{mKrig.coef}
\alias{mKrig.trace}
\alias{print.mKrig}
\alias{print.mKrigSummary}
\alias{summary.mKrig}
\alias{mKrigCheckXY}

\title{"micro Krig"  Spatial process estimate of a curve or surface, 
"kriging" with a known covariance function. }
\description{
This is a simple version of the Krig function that is 
optimized for large data sets, sparse linear algebra, and a clear exposition of the
computations. Lambda, the smoothing parameter must be fixed. 
This function is called higher level functions for maximum likelihood estimates of 
covariance paramters. 
}
\usage{
mKrig(x, y, weights = rep(1, nrow(x)), Z = NULL,
                 cov.function = "stationary.cov", cov.args = NULL,
                 lambda = 0, m = 2, chol.args = NULL, find.trA = TRUE,
                 NtrA = 20, iseed = 123, llambda = NULL, na.rm = FALSE,
                  collapseFixedEffect = TRUE,
                 ...)
  
\method{predict}{mKrig}( object, xnew=NULL,ynew=NULL, grid.list = NULL,
derivative=0,
Z=NULL,drop.Z=FALSE,just.fixed=FALSE,
collapseFixedEffect = object$collapseFixedEffect, ...)

\method{summary}{mKrig}(object, ...)

\method{print}{mKrig}( x, digits=4,... )
\method{print}{mKrigSummary}( x, digits=4,... )

mKrig.coef(object, y, collapseFixedEffect=TRUE)

mKrig.trace( object, iseed, NtrA)

mKrigCheckXY(x, y, weights, Z, na.rm) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{collapseFixedEffect}{ If replicated fields are given to mKrig (i.e.
\code{y} has more than one column) there is the choice of estimating the
fixed effect coefficients (\code{d} in the returned object) separately
for each replicate or pooling across replicates and deriving a single
estimate. If \code{collapseFixedEffect} is TRUE (default) the estimates are
pooled. }
  
\item{chol.args}{A list of optional arguments (pivot, nnzR) that will
be used with the call to the cholesky decomposition. Pivoting is done
by default to make use of sparse matrices when they are
generated. This argument is useful in some cases for sparse covariance
functions to reset the memory parameter nnzR.  (See example below.)}

\item{cov.args}{A list of optional arguments that will be used in
  calls to the covariance function.}

\item{cov.function}{The name, a text string of the covariance function.}

\item{derivative}{If zero the surface will be evaluated. If not zero
  the matrix of partial derivatives will be computed.}

\item{digits}{Number of significant digits used in printed output.}

\item{drop.Z}{If true the fixed part will only be evaluated at the
polynomial part of the fixed model. The contribution from the other
covariates will be omitted.}

\item{find.trA}{
If TRUE will estimate the effective degrees of freedom using 
  a simple Monte Carlo method. This will add to the computational 
  burden by approximately NtrA solutions of the linear system but 
  the cholesky decomposition is reused.}

\item{grid.list}{A grid.list to evaluate the surface in place of specifying 
  arbitrary locations.}

\item{iseed}{Random seed ( using \code{set.seed(iseed)}) used to generate 
  iid normals for Monte Carlo estimate of the trace.}

\item{just.fixed}{If TRUE only the predictions for the fixed part of
  the model will be evaluted.}
  
\item{lambda}{ Smoothing parameter or equivalently the ratio between 
  the nugget and process varainces.}

\item{llambda}{If not \code{NULL} then \code{lambda = exp( llambda)}}

\item{m}{ The degree of the polynomial used in teh fixed part is (m-1)}

\item{na.rm}{If TRUE NAs in y are omitted along with corresonding rows of x.}

\item{NtrA}{Number of Monte Carlo samples for the trace. But if NtrA is 
  greater than or equal to the number of observations the trace 
  is computed exactly.}

\item{object}{Object returned by mKrig. (Same as "x"
  in the print function.)}
  
\item{weights}{Precision ( 1/variance) of each observation}

\item{x}{Matrix of unique spatial locations (or in print or surface 
the returned mKrig object.)}

\item{xnew}{Locations for predictions.}

\item{y}{ Vector or matrix of observations at spatial locations, 
  missing values are not allowed! Or in mKrig.coef a new 
  vector of observations. If y is a matrix the columns are 
  assumed to be independent replicates of the spatial field.  I.e.
  observation vectors generated 
  from the same covariance and measurment error model but 
  independent from each other. 
}

\item{ynew}{New observation vector. \code{mKrig} will reuse matrix
  decompositions and find the new fit to these data.}

\item{Z}{ Linear covariates to be included in fixed part of the 
  model that are distinct from the default low order 
  polynomial in \code{x}. (NOTE the order of the polynomial
  determined by \code{m})}

\item{\dots}{ In \code{mKrig} and \code{predict} additional arguments
that will be passed to the covariance function.}

}
\details{
 This function is an abridged version of Krig. The m stand for micro
and this function focuses on the computations in Krig.engine.fixed
done for a fixed lambda parameter, for unique spatial locations and
for data without missing values.

These restrictions simplify the code for reading. Note that also
little checking is done and the spatial locations are not transformed
before the estimation.  Because most of the operations are linear
algebra this code has been written to handle multiple data
sets. Specifically if the spatial model is the same except for
different observed values (the y's), one can pass \code{y} as a matrix
and the computations are done efficiently for each set.  Note that
this is not a multivariate spatial model just an efficient computation
over several data vectors without explicit looping.A big difference in
the computations is that an exact expression for thetrace of the
smoothing matrix is (trace A(lambda)) is computationally expensive and
a Monte Carlo approximation is supplied instead.

See \code{predictSE.mKrig} for prediction standard errors and 
\code{sim.mKrig.approx} to quantify the uncertainty in the estimated function using conditional 
simulation.  

\code{predict.mKrig} will evaluate the derivatives of the estimated
function if derivatives are supported in the covariance function.  For
example the wendland.cov function supports derivatives.


\code{summary.mKrig} creates a list of class \code{mKrigSummary} along with a table of standard errors for the fixed linear parameters.

\code{print.mKrigSummary} prints the  \code{mKrigSummary} object and adds some more explanation about the model and results

\code{print.mKrig} prints a summary for the \code{mKrig} object that the combines the summary and  print methods.


\code{mKrig.coef} finds the "d" and "c" coefficients represent the
solution using the previous cholesky decomposition for a new data
vector. This is used in computing the prediction standard error in
predictSE.mKrig and can also be used to evalute the estimate
efficiently at new vectors of observations provided the locations and
covariance remain fixed.

Sparse matrix methods are handled through overloading the usual linear
algebra functions with sparse versions. But to take advantage of some
additional options in the sparse methods the list argument chol.args
is a device for changing some default values. The most important of
these is \code{nnzR}, the number of nonzero elements anticipated in
the Cholesky factorization of the postive definite linear system used
to solve for the basis coefficients. The sparse of this system is
essentially the same as the covariance matrix evalauted at the
observed locations.  As an example of resetting \code{nzR} to 450000
one would use the following argument for chol.args in mKrig:

\code{ chol.args=list(pivot=TRUE,memory=list(nnzR= 450000))}

\code{mKrig.trace} This is an internal function called by \code{mKrig}
to estimate the effective degrees of freedom.  The Kriging surface
estimate at the data locations is a linear function of the data and
can be represented as A(lambda)y.  The trace of A is one useful
measure of the effective degrees of freedom used in the surface
representation. In particular this figures into the GCV estimate of
the smoothing parameter.  It is computationally intensive to find the
trace explicitly but there is a simple Monte Carlo estimate that is
often very useful.  If E is a vector of iid N(0,1) random variables
then the trace of A is the expected value of t(E)AE. Note that AE is
simply predicting a surface at the data location using the synthetic
observation vector E. This is done for \code{NtrA} independent N(0,1)
vectors and the mean and standard deviation are reported in the
\code{mKrig} summary.  Typically as the number of observations is
increased this estimate becomse more accurate. If NtrA is as large as
the number of observations (\code{np}) then the algorithm switches to
finding the trace exactly based on applying A to \code{np} unit
vectors.
}
\value{
  \item{d}{Coefficients of the polynomial fixed part and if present
  the covariates (Z).To determine which is which the logical vector
  ind.drift also part of this object is TRUE for the polynomial
  part. }

  \item{c}{ Coefficients of the nonparametric part.}

  \item{nt}{ Dimension of fixed part.}

  \item{np}{ Dimension of c.}

  \item{nZ}{Number of columns of Z covariate matrix (can be zero).}

   \item{ind.drift}{Logical vector that indicates polynomial
   coefficients in the \code{d} coefficients vector. This is helpful
   to distguish between polynomial part and the extra covariates
   coefficients associated with Z. }
 
 \item{lambda.fixed}{The fixed lambda value}

  \item{x}{Spatial locations used for fitting.}

  \item{knots}{The same as x}

  \item{cov.function.name}{Name of covariance function used.}

 \item{args}{ A list with all the covariance arguments that were
 specified in the call.}

 \item{m}{Order of fixed part polynomial.}

  \item{chol.args}{ A list with all the cholesky arguments that were
         specified in the call.}

 \item{call}{ A copy of the call to mKrig.}

 \item{non.zero.entries}{ Number of nonzero entries in the covariance
matrix for the process at the observation locations.}

 \item{shat.MLE}{MLE of sigma.}

 \item{rho.MLE}{MLE or rho.}

 \item{rhohat}{Estimate for rho adjusted for fixed model degrees of
 freedom (ala REML).}

 \item{lnProfileLike}{log Profile likelihood for lambda}

 \item{lnDetCov}{Log determinant of the covariance matrix for the
  observations having factored out rho.}

 \item{Omega}{GLS covariance for the estimated parameters in the fixed
part of the model (d coefficients0.}

 \item{qr.VT, Mc}{QR and cholesky matrix decompositions needed to
recompute the estimate for new observation vectors.}

 \item{fitted.values, residuals}{Usual predictions from fit.}

 \item{eff.df}{Estimate of effective degrees of freedom. Either the
mean of the Monte Carlo sample or the exact value. }

 \item{trA.info}{If NtrA ids less than \code{np} then the individual
members of the Monte Carlo sample and \code{sd(trA.info)/ sqrt(NtrA)}
is an estimate of the standard error. If NtrA is greater than or equal
to \code{np} then these are the diagonal elements of A(lamdba).}

 \item{GCV}{Estimated value of the GCV function.}

 \item{GCV.info}{Monte Carlo sample of GCV functions}

}
\author{Doug Nychka, Reinhard Furrer, John Paige}
\seealso{ Krig, surface.mKrig, Tps, fastTps, predictSurface, predictSE.mKrig, sim.mKrig.approx, 
  \code{ \link{mKrig.grid}}}
\examples{
#
# Midwest ozone data  'day 16' stripped of missings 
  data( ozone2)
  y<- ozone2$y[16,]
  good<- !is.na( y)
  y<-y[good]
  x<- ozone2$lon.lat[good,]
# nearly interpolate using defaults (Exponential covariance range = 2.0)
# see also mKrigMLEGrid to choose lambda by maxmimum likelihood
  out<- mKrig( x,y, theta = 2.0, lambda=.01)
  out.p<- predictSurface( out)
  surface( out.p)
#
# NOTE this should be identical to 
# Krig( x,y, theta=2.0, lambda=.01) 

##############################################################################
# an example using a "Z" covariate and the Matern family
#  again see mKrigMLEGrid to choose parameters by MLE.
data(COmonthlyMet)
yCO<- CO.tmin.MAM.climate
good<- !is.na( yCO)
yCO<-yCO[good]
xCO<- CO.loc[good,]
Z<- CO.elev[good]
out<- mKrig( xCO,yCO, Z=Z, cov.function="stationary.cov", Covariance="Matern",
                    theta=4.0, smoothness=1.0, lambda=.1)
set.panel(2,1)
# quilt.plot with elevations
quilt.plot( xCO, predict(out))
# Smooth surface without elevation linear term included
surface( out)
set.panel()

#########################################################################
# Interpolate using tapered version of the exponential, 
# the taper scale is set to 1.5 default taper covariance is the Wendland.
# Tapering will done at a scale of 1.5 relative to the scaling 
# done through the theta  passed to the covariance function.
data( ozone2)
  y<- ozone2$y[16,]
  good<- !is.na( y)
  y<-y[good]
  x<- ozone2$lon.lat[good,]
  mKrig( x,y,cov.function="stationary.taper.cov",
       theta = 2.0, lambda=.01, 
       Taper="Wendland",  Taper.args=list(theta = 1.5, k=2, dimension=2)
           ) -> out2

# Try out GCV on a grid of lambda's.
# For this small data set 
# one should really just use Krig or Tps but this is an example of
# approximate GCV that will work for much larger data sets using sparse 
# covariances and the Monte Carlo trace estimate
#
# a grid of lambdas:
  lgrid<- 10**seq(-1,1,,15) 
  GCV<- matrix( NA, 15,20)
  trA<-  matrix( NA, 15,20)
  GCV.est<- rep( NA, 15)
  eff.df<- rep( NA, 15)
  logPL<- rep( NA, 15) 
# loop over lambda's
  for(  k in 1:15){
      out<- mKrig( x,y,cov.function="stationary.taper.cov",
                    theta = 2.0, lambda=lgrid[k],
          Taper="Wendland",  Taper.args=list(theta = 1.5, k=2, dimension=2)  ) 
      GCV[k,]<- out$GCV.info
      trA[k,]<- out$trA.info
      eff.df[k]<- out$eff.df
      GCV.est[k]<- out$GCV
      logPL[k]<- out$lnProfileLike
  }
#
# plot the results different curves are for individual estimates  
# the two lines are whether one averages first the traces or the GCV criterion.
#
  par( mar=c(5,4,4,6))
  matplot( trA, GCV, type="l", col=1, lty=2,
            xlab="effective degrees of freedom", ylab="GCV")
  lines( eff.df, GCV.est, lwd=2, col=2)
  lines( eff.df, rowMeans(GCV), lwd=2)
# add exact GCV computed by Krig 
  out0<-  Krig( x,y,cov.function="stationary.taper.cov",
              theta = 2.0, 
              Taper="Wendland", Taper.args=list(theta = 1.5, k=2, dimension=2),
              spam.format=FALSE)  
  lines( out0$gcv.grid[,2:3], lwd=4, col="darkgreen")

# add profile likelihood 
  utemp<- par()$usr
  utemp[3:4] <- range( -logPL)
  par( usr=utemp)
  lines( eff.df, -logPL, lwd=2, col="blue", lty=2)
  axis( 4)
  mtext( side=4,line=3, "-ln profile likelihood", col="blue")
  title( "GCV ( green = exact) and  -ln profile likelihood", cex=2)

#########################################################################
# here is a series of examples with bigger datasets  
# using a compactly supported covariance directly

set.seed( 334)
N<- 1000
x<- matrix( 2*(runif(2*N)-.5),ncol=2)
y<- sin( 1.8*pi*x[,1])*sin( 2.5*pi*x[,2]) + rnorm( 1000)*.1
  
look2<-mKrig( x,y, cov.function="wendland.cov",k=2, theta=.2, 
            lambda=.1)

# take a look at fitted surface
predictSurface(look2)-> out.p
surface( out.p)

# this works because the number of nonzero elements within distance theta
# are less than the default maximum allocated size of the 
# sparse covariance matrix. 
#  see  options() for the default values. The names follow the convention
# spam.arg where arg is the name of the spam component 
#   e.g. spam.nearestdistnnz

# The following will give a warning for theta=.9 because 
# allocation for the  covariance matirx storage is too small. 
# Here theta controls the support of the covariance and so 
# indirectly the  number of nonzero elements in the sparse matrix

\dontrun{
 look2<- mKrig( x,y, cov.function="wendland.cov",k=2, theta=.9, lambda=.1)
}

# The warning resets the memory allocation  for the covariance matrix
# according the to values   options(spam.nearestdistnnz=c(416052,400))'
# this is inefficient becuase the preliminary pass failed. 

# the following call completes the computation in "one pass"
# without a warning and without having to reallocate more memory. 

options( spam.nearestdistnnz=c(416052,400))
  look2<- mKrig( x,y, cov.function="wendland.cov",k=2,
                    theta=.9, lambda=1e-2)
# as a check notice that 
#   print( look2)
# reports the number of nonzero elements consistent with the specifc allocation
# increase in spam.options


# new data set of 1500 locations
  set.seed( 234)
  N<- 1500
  x<- matrix( 2*(runif(2*N)-.5),ncol=2)
  y<- sin( 1.8*pi*x[,1])*sin( 2.5*pi*x[,2]) + rnorm( N)*.01

\dontrun{  
# the following is an example of where the allocation  (for nnzR) 
# for the cholesky factor is too small. A warning is issued and 
# the allocation is increased by 25% in this example
#
 look2<- mKrig( x,y, 
            cov.function="wendland.cov",k=2, theta=.1, lambda=1e2  )
}
# to avoid the warning 
 look2<-mKrig( x,y, 
            cov.function="wendland.cov", k=2, theta=.1, 
            lambda=1e2, chol.args=list(pivot=TRUE, memory=list(nnzR= 450000)))

###############################################################################
# fiting multiple data sets
#
#\dontrun{ 
  y1<- sin( 1.8*pi*x[,1])*sin( 2.5*pi*x[,2]) + rnorm( N)*.01
  y2<- sin( 1.8*pi*x[,1])*sin( 2.5*pi*x[,2]) + rnorm( N)*.01
  Y<- cbind(y1,y2)
  look3<- mKrig( x,Y,cov.function="wendland.cov",k=2, theta=.1, 
            lambda=1e2  )
# note slight difference in summary because two data sets have been fit.
  print( look3)
#}

##################################################################
# finding a good choice for theta as a taper 

# Suppose the target is a spatial prediction using roughly 50 nearest neighbors
# (tapering covariances is effective for roughly 20 or more in the situation of 
#  interpolation) see Furrer, Genton and Nychka (2006).
# take a look at a random set of 100 points to get idea of scale
# and saving  computation time by not  looking at the complete set
# of points
# NOTE: This could also be done directly using the  FNN package for finding nearest 
# neighbors
  set.seed(223)
  ind<- sample( 1:N,100)
  hold<- rdist( x[ind,], x)
  dd<- apply( hold, 1, quantile, p= 50/N )
  dguess<- max(dd)
# dguess is now a reasonable guess at finding cutoff distance for
# 50 or so neighbors
# full distance matrix excluding distances greater than dguess
  hold2<- nearest.dist( x, x, delta= dguess )
# here is trick to find the number of nonsero rows for a matrix in spam format. 
  hold3<-  diff( hold2@rowpointers)
#  min( hold3) = 43   which we declare close enough. This also counts the diagonal
# So there are a minimum of 42 nearest neighbors  ( median is 136)
# see  table( hold3) for the distribution 
# now the following will use no less than 43 - 1  nearest neighbors 
# due to the tapering. 
\dontrun{
  mKrig( x,y, cov.function="wendland.cov",k=2, theta=dguess, 
            lambda=1e2) ->  look2
}

###############################################################################
# use precomputed distance matrix
#
\dontrun{ 
  y1<- sin( 1.8*pi*x[,1])*sin( 2.5*pi*x[,2]) + rnorm( N)*.01
  y2<- sin( 1.8*pi*x[,1])*sin( 2.5*pi*x[,2]) + rnorm( N)*.01
  Y<- cbind(y1,y2)
  #precompute distance matrix in compact form
  distMat = rdist(x, compact=TRUE)
  look3<- mKrig( x,Y,cov.function="stationary.cov", theta=.1, 
            lambda=1e2, distMat=distMat )
  #precompute distance matrix in standard form
  distMat = rdist(x)
  look3<- mKrig( x,Y,cov.function="stationary.cov", theta=.1, 
            lambda=1e2, distMat=distMat )
}
}
\references{
\url{https://github.com/NCAR/Fields}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial }
