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

\name{qsreg}
\alias{qsreg}
\title{
  Quantile or Robust spline regression  
}
\description{
Uses a penalized likelihood approach to estimate the conditional  
quantile function for regression data. This method is only implemented  
for univariate data. For the pairs (X,Y) the  
conditional quantile, f(x), is  P( Y<f(x)| X=x) = alpha. This estimate  
is useful for determining the envelope of a scatterplot or assessing  
departures from a constant variance with respect to the independent  
variable.    
}
\usage{
qsreg(x, y, lam = NA, maxit = 50, maxit.cv = 10, tol =
                 1e-07, offset = 0, sc = sqrt(var(y)) * 1e-05, alpha =
                 0.5, wt = rep(1, length(x)), cost = 1, nstep.cv = 80,
                 hmin = NA, hmax = NA, trmin = 2 * 1.05, trmax = 0.95
                 * length(unique(x)))
}
\arguments{
\item{x}{
Vector of the independent variable in  y = f(x) + e}
\item{y}{
Vector of the dependent variable}
\item{lam}{
Values of the smoothing parameter. If omitted is found by GCV based on the  
the quantile criterion 
}
\item{maxit}{
Maximum number of iterations used to estimate each quantile spline. 
}
\item{maxit.cv}{
Maximum number of iterations to find GCV minimum. 
}
\item{tol}{
Tolerance for convergence when computing quantile spline. 
}
\item{cost}{
Cost value used in the GCV criterion. Cost=1 is the usual GCV  
denominator. 
}
\item{offset}{
Constant added to the effective degrees of freedom in the GCV function.  
}
\item{sc}{
Scale factor for rounding out the absolute value function at zero to a 
quadratic. Default is a small scale to produce something more like 
quantiles. Scales on the order of the residuals will result is a robust 
regression fit using the Huber weight function. The default is 1e-5 of the 
variance of the Y's. The larger this value the better behaved the problem 
is numerically and requires fewer iterations for convergence at each new 
value of lambda. 
}
\item{alpha}{
Quantile to be estimated. Default is find the median. 
}
\item{wt}{
Weight vector default is constant values. Passing nonconstant weights is a 
pretty strange thing to do.  
}
\item{nstep.cv}{
Number of points used in CV grid search 
}
\item{hmin}{
Minimum value of log( lambda) used for GCV grid search. 
}
\item{hmax}{
Maximum value of log( lambda) used for GCV grid search. 
}
\item{trmin}{
Minimum value of effective degrees of freedom in model used 
for specifying the range of lambda in the GCV grid search. 
}
\item{trmax}{
Maximum value of effective degrees of freedom in model used 
for specifying the range of lambda in the GCV grid search. 
}

}
\value{

\item{trmin trmax }{
Define the minimum and maximum values for the CV grid search in terms of 
the effective number of parameters. (see hmin, hmax) 
Object of class qsreg with many arguments similar to a sreg object.  
One difference is that cv.grid has five columns the last being  
the number of iterations for convergence at each value of lambda.  
}
}
\details{
This is an experimental function to find the smoothing parameter for a  
quantile or robust spline using a more appropriate criterion than mean squared  
error prediction.  
The quantile spline is found by an iterative algorithm using weighted  
least squares cubic splines. At convergence the estimate will also be a  
weighted natural  cubic spline but the weights will depend on the 
estimate. 
Alternatively at convergence the estimate will be a least squares spline applied to the 
empirical psuedo data. The user is referred to the paper by Oh and Nychka ( 2002) for the 
details and properties of the robust cross-validation using empirical psuedo data.
Of course these weights are crafted so that the resulting spline is an  
estimate of the alpha quantile instead of the mean. CV as function of 
lambda can be strange so it should be plotted. 
  
}
\seealso{
\code{\link{sreg}} 
}
\examples{

     # fit a CV  quantile spline
     fit50<- qsreg(rat.diet$t,rat.diet$con)
     # (default is .5 so this is an estimate of the conditional median)
     # control group of rats.
     plot( fit50)
     predict( fit50)
     # predicted values at data points
     xg<- seq(0,110,,50)
     plot( fit50$x, fit50$y)
     lines( xg, predict( fit50, xg))

     # A robust fit to rat diet data
     # 
     SC<- .5* median(abs((rat.diet$con- median(rat.diet$con))))
     fit.robust<- qsreg(rat.diet$t,rat.diet$con, sc= SC)
     plot( fit.robust)

     # The global GCV function suggests little smoothing so 
     # try the local
     # minima with largest lambda instead of this default value.
     # one should should consider redoing the three quantile fits in this
     # example after looking at the cv functions and choosing a good value for
     #lambda
     # for example
     lam<- fit50$cv.grid[,1]
     tr<- fit50$cv.grid[,2]
     # lambda close to df=6
     lambda.good<- max(lam[tr>=6])
     fit50.subjective<-qsreg(rat.diet$t,rat.diet$con, lam= lambda.good)
     fit10<-qsreg(rat.diet$t,rat.diet$con, alpha=.1, nstep.cv=200)
     fit90<-qsreg(rat.diet$t,rat.diet$con, alpha=.9, nstep.cv=200)
     # spline fits at 50 equally spaced points
     sm<- cbind(
 
     predict( fit10, xg),
     predict( fit50.subjective, xg),predict( fit50, xg),
     predict( fit90, xg))
 
     # and now zee data ...
     plot( fit50$x, fit50$y)
     # and now zee quantile splines at 10% 50% and 90%.
     #
     matlines( xg, sm, col=c( 3,3,2,3), lty=1) # the spline
  
}
\keyword{smooth}
% docclass is function
% Converted by Sd2Rd version 1.21.
