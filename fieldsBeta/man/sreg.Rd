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

\name{sreg}
\alias{sreg}
\alias{predict.sreg}
\title{
  Cubic smoothing spline regression
}
\description{
Fits a cubic smoothing spline to univariate data. The amount of
smoothness can be specified or estimated from the data by GCV. 
<!--brief description-->
}
\usage{
sreg(x, y, lambda = NA, df = NA, offset = 0, weights =
                 rep(1, length(x)), cost = 1, nstep.cv = 80, tol =
                 1e-05, find.diagA = TRUE, trmin = 2.01, trmax = NA,
                 lammin = NA, lammax = NA, verbose = FALSE, do.cv =
                 TRUE, method = "GCV", rmse = NA, na.rm = TRUE, digits
                 = 8)


\method{predict}{sreg}(object, x, derivative = 0, model = 1,...) 
}

\arguments{
\item{x}{
Vector of x value}

\item{y}{
Vector of y values}

\item{lambda}{
Single smoothing parameter or a vector of values . If omitted 
smoothing parameter estimated by GCV. NOTE: lam here is equivalent to 
the value lambda*N  in Tps/Krig  where N is the number of unique observations.
See example below.}

\item{object}{An sreg object.}
\item{derivative}{Order of deriviatve to evaluate. Must be 0,1, or 2.}
\item{df}{
Amount of smoothing in term of effective degrees of freedom for the
spline}

\item{offset}{
an offset added to the term cost*degrees of freedom in the denominator of
the GCV function. (This would be used for adjusting the df from fitting
other models such as in back-fitting additive models.)}

\item{model}{Specifies which model parameters to use.}

\item{weights}{
A vector that is proportional to the reciprocal variances of the
errors.}

\item{cost}{
Cost value to be used in the GCV criterion.}

\item{nstep.cv }{
Number of grid points of smoothing parameter for GCV grid search.}

\item{tol}{Tolerance for convergence in minimizing the GCV or other
criteria to estimate the smoothing parameter.}
\item{find.diagA}{
If TRUE  calculates the diagonal elements of the smoothing matrix. The
effective
number of degrees of freedom is the sum of these diagonal elements.
Default is true. This requires more stores if a grid of smoothing
parameters is passed. ( See returned values below.)}

\item{trmin}{
Sets the minimum of the smoothing parameter range  for the GCV grid
search in terms of effective degrees of freedom.}

\item{trmax}{
Sets the maximum of the smoothing parameter range  for the GCV grid
search in terms of effective degrees of freedom. If NA the range is set
to .99 of number of unique locations.}

\item{lammin}{
Same function as trmin but in the lambda scale.}

\item{lammax}{
Same function as trmax but in the lambda scale.}

\item{verbose}{
Print out all sorts of debugging info. Default is falseof course!}

\item{do.cv}{
Evaluate the spline at the GCV minimum. Default is true.}

\item{method}{
A character string giving the 
method for determining the smoothing
parameter. Choices are
"GCV", "GCV.one", "GCV.model", "pure error", "RMSE". Default is "GCV". }

\item{rmse}{
Value of the root mean square error to match by varying lambda.}

\item{na.rm}{If TRUE NA's are removed from y before analysis.}

\item{digits}{Number of significant digits used to determine replicate x values.}

\item{\dots}{Other optional arguments to pass to the predict function.}

}
\value{
Returns a list of class sreg. 
Some of the returned components are 

\item{call}{
Call to the function 
}
\item{yM}{
Vector of dependent variables. If replicated data is given these are the
replicate group means. }

\item{xM}{
Unique x values matching the y's. }
\item{weights}{
Proportional to reciprocal variance of each data point.
}
\item{weightsM}{
Proportional to reciprocal pooled variance of each 
replicated mean data value (xM).}

\item{x}{
Original x data. }

\item{y}{
Original y data. }

\item{method}{
Method used to find the smoothing parameter.}

\item{pure.ss}{
Pure error sum of squares from replicate groups. }

\item{tauHat.pure.error}{
Estimate of tau from replicate groups.}

\item{tauHat.GCV}{
Estimate of tau using estimated lambda from GCV minimization }

\item{trace}{
Effective degrees of freedom for the spline estimate(s)}

\item{gcv.grid}{
Values of trace, GCV, tauHat. etc. for a grid of smoothing parameters.
If lambda ( or df) is specified those values are used. }

\item{lambda.est}{
Summary of various estimates of the smoothing parameter}

\item{lambda}{
If lambda is specified the passed  vector, if missing the estimated value.}

\item{residuals}{
Residuals from spline(s). If lambda or df is specified the residuals from
these values. If lambda and df are omitted then the spline having
estimated lambda. This will be a matrix with as many columns as the values
of lambda. }

\item{fitted.values}{
Matrix of fitted values. See notes on residuals. }

\item{predicted}{
A list with components  x and y. x is the unique values of xraw in sorted
order. y is a matrix of the spline estimates at these values. }

\item{eff.df}{
Same as trace.}

\item{diagA}{
Matrix containing diagonal elements of the smoothing matrix. Number of
columns is the number of lambda values. 
WARNING: If there is replicated data the
diagonal elements are those for the smoothing the group means at the
unique x locations. }
}

\details{
MODEL: The assumed model is Y.k=f(x.k) +e.k where e.k should be
approximately
normal and independent errors with variances tau**2/w.k

ESTIMATE: A smoothing spline is a locally weighted average of the y's
based 
on the relative locations of the x values. Formally the estimate is 
the curve that minimizes the criterion: 

    
(1/n) sum(k=1,n) w.k( Y.k - f( X.k))**2  + lambda R(f) 

where R(f) is the integral of the squared second derivative of f over 
the range of the X values. Because of the inclusion of the (1/n) in the
sum of squares the lambda parameter in sreg corresponds to the a value of
lambda*n in the Tps function and in the Krig function. 
 
 The solution to this minimization is a piecewise cubic polynomial with
the join points at the unique set of X values. The polynomial segments
are constructed so that the entire curve has continuous first and second
derivatives and the second and third derivatives are zero at the
boundaries.  The smoothing has the range [0,infinity]. Lambda equal to
zero gives a cubic spline interpolation of the data. As lambda diverges
to infinity ( e.g lambda =1e20) the estimate will converge to the
straight line estimated by least squares.

    The values of the estimated function at the data points can be
expressed in the matrix form:
  
    predicted values= A(lambda)Y 

where A is an nXn symmetric matrix that does NOT depend on Y. 
The diagonal elements are the leverage values for the estimate and the 
sum of these  (trace(A(lambda)) can be interpreted as the effective 
number of parameters that are used to define the spline function. 
IF there are replicate points the A matrix is the result of finding group
averages and applying a weighted spline to the means. 
The A matrix is also used to find "Bayesian" confidence intervals for the 
estimate, see the example below. 

CROSS-VALIDATION:The GCV criterion with no replicate points for a fixed
value of lambda is

 (1/n)(Residual sum of squares)/((1-(tr(A)-offset)*cost + offset)/n)**2, 

Usually offset =0 and cost =1. Variations on GCV with replicate points are
described in the documentation help file for Krig.  With an appropriate
choice for the smoothing parameter, the estimate of tau**2 is found by
(Residual sum of squares)/tr(A).

COMPUTATIONS: The computations for 1-d splines exploit the banded
structure of the matrices needed to solve for the spline coefficients.
Banded structure also makes it possible to get the diagonal elements of A
quickly. This approach is different from the algorithms in Tps and
tremendously more efficient for larger numbers of unique x values ( say >
200). The advantage of Tps is getting "Bayesian" standard errors at
predictions different from the observed x values. This function is similar
to the S-Plus smooth.spline. The main advantages are more information and
control over the choice of lambda and also the FORTRAN source code is
available (css.f). 

See also the function \code{splint} which is designed to be a bare bones
but fast smoothing spline. 
 
}
\seealso{
Krig, Tps, splint
}
\examples{
# fit a GCV spline to  
# control group of rats.  
fit<- sreg(rat.diet$t,rat.diet$con)
summary( fit)

set.panel(2,2)
plot(fit)                       # four diagnostic plots of  fit 
set.panel()

predict( fit) # predicted values at data points 

xg<- seq(0,110,,50) 
sm<-predict( fit, xg) # spline fit at 50 equally spaced points 
der.sm<- predict( fit, xg, deriv=1) # derivative of spline fit 
set.panel( 2,1) 
plot( fit$x, fit$y) # the data 
lines( xg, sm) # the spline 
plot( xg,der.sm, type="l") # plot of estimated derivative 
set.panel() # reset panel to 1 plot


# the same fit using  the thin plate spline numerical algorithms 
# sreg does not scale the obs so instruct Tps not to sacel either
# this will make lambda comparable within factor of n. 

   fit.tps<-Tps( rat.diet$t,rat.diet$con, scale="unscaled")
   summary( fit.tps) 

# compare sreg and Tps results to show the adjustment to lambda.

   predict( fit)-> look
   predict( fit.tps, lambda=fit$lambda*fit$N)-> look2
   test.for.zero( look, look2) # silence means it checks to 1e-8

# finding approximate standard errors at observations

SE<- fit$tauHat.GCV*sqrt(fit$diagA)

# compare to predictSE( fit.tps) differences are due to 
# slightly different lambda values and using tauHat.MLE instad of tauHat.GCV
#

# 95% pointwise prediction intervals
Zvalue<-  qnorm(.0975)
upper<- fit$fitted.values + Zvalue* SE
lower<- fit$fitted.values - Zvalue* SE
#
# conservative, simultaneous Bonferroni bounds
#
ZBvalue<-  qnorm(1- .025/fit$N)
upperB<- fit$fitted.values + ZBvalue* SE
lowerB<- fit$fitted.values - ZBvalue* SE
#
# take a look

plot( fit$x, fit$y)
lines( fit$predicted, lwd=2)
matlines( fit$x, 
cbind( lower, upper, lowerB, upperB), type="l", col=c( 2,2,4,4), lty=1)
title( "95 pct pointwise  and simultaneous intervals")
# or try the more visually  honest:
plot( fit$x, fit$y)
lines( fit$predicted, lwd=2)
segments(  fit$x, lowerB, fit$x, upperB, col=4)
segments(  fit$x, lower, fit$x, upper, col=2, lwd=2)
title( "95 pct pointwise  and simultaneous intervals")

set.panel( 1,1)
}
\keyword{smooth}
