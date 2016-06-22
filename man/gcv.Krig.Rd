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

\name{gcv.Krig}
\alias{gcv.Krig}
\alias{gcv.sreg}

\title{Finds profile likelihood and GCV estimates of 
smoothing parameters for splines and Kriging.}
\description{
 This is a secondary function that will use the computed Krig object and
find various estimates of the smoothing parameter lambda. These are
several different flavors of cross-validation, a moment matching
strategy and the profile likelihood.  This function can also be used
independently with different data sets (the y's) if the covariates ( the
x's) are the same and thus reduce the computation. 
}
\usage{
gcv.Krig(
out, lambda.grid = NA, cost = 1, nstep.cv = 200, rmse
                 = NA, verbose = FALSE, tol = 1e-05, offset = 0, y =
                 NULL, give.warnings = TRUE)

gcv.sreg (
out, lambda.grid = NA, cost = 1, nstep.cv = 80, rmse =
                 NA, offset = 0, trmin = NA, trmax = NA, verbose =
                 FALSE, tol = 1e-05, give.warnings = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out}{ A Krig or sreg object.}
  \item{lambda.grid}{ Grid of lambdas for coarse search. The default is 
equally spaced on effective degree of freedom scale. }
  \item{cost}{ Cost used in GCV denominator }
  \item{nstep.cv}{ Number of grid points in coarse search. }
  \item{rmse}{ Target root mean squared error to match with 
                  the estimate of sigma**2 }
  \item{verbose}{ If true prints intermediate results.  }
  \item{tol}{ Tolerance in delcaring convergence of golden section search or bisection search. }
  \item{offset}{ Additional degrees of freedom to be added into the GCV denominator.
}
  \item{y}{A new data vector to be used in place of the one associated with the
Krig object (obj) }
  \item{give.warnings}{ If FALSE will suppress warnings about grid search being out of 
range for various estimates based on GCV and REML.}
 
\item{trmin}{Minimum value of lambda for grid search specified in terms
of effective degrees of freedom.}
\item{trmax}{Maximum value for grid search.}

}
\details{

 This function finds several estimates of the smoothing parameter using
first a coarse grid search followed by a refinement using a minimization (
in the case of GCV or maximum likelihood) or bisection in the case of
mathcing the rmse. Details of the estimators can be found in the help file
for the Krig function. 

The Krig object passed to this function has some matrix decompostions that
facilitate rapid computation of the GCV and ML functions and do not depend
on the independent variable. This makes it possible to compute the Krig
object once and to reuse the decompostions for multiple data sets. (But
keep in mind if the x values change then the object must be recalculated.)
The example below show show this can be used for a simulation study on the
variability for estimating the smoothing parameter. 


}
\value{A list giving a summary of estimates and diagonostic details with the 
following components:
  \item{gcv.grid }{ A matrix describing results of the 
coarse search rows are values of lambda and the columns are 
lambda= value of smoothing parameter, 
trA=effective degrees of freedom, 
GCV=Usual GCV criterion, 
GCV.one=GCV criterion leave-one-out, 
GCV.model= GCV based on average response in the case of replicates, 
shat= Implied estimate of sigma ,
-Log Profile= negative log of profiel likelihood for the lambda. 
}
  \item{lambda.est}{Summary table of all estimates
Rows index different types of estimates: 
GCV, GCV.model, GCV.one, RMSE, pure error, -Log Profile
and the columns are the estimated values  for lambda, trA, GCV, shat.
}
}
\author{Doug Nychka}

\seealso{ 
\code{\link{Krig}},
\code{\link{Tps}},
\code{\link{predict.Krig}} }
\examples{

# 
Tps( ChicagoO3$x, ChicagoO3$y)-> obj # default is to find lambda by GCV
summary( obj)

gcv.Krig( obj)-> out
print( out$lambda.est) # results agree with Tps summary

sreg( rat.diet$t, rat.diet$trt)-> out
gcv.sreg( out, tol=1e-10) # higher tolerance search for minimum 
\dontrun{
# a simulation example
x<- seq( 0,1,,150)
f<-  x**2*( 1-x)
f<- f/sqrt( var( f))

set.seed(123) # let's all use the same seed
sigma<- .1
y<- f + rnorm( 150)*sigma

Tps( x,y)-> obj # create Krig object

hold<- hold2<- matrix( NA, ncol=6, nrow=200)

for( k in 1:200){
# look at GCV estimates of lambda
# new data simulated
   y<- f + rnorm(150)*sigma 
# save GCV estimates
lambdaTable<- gcv.Krig(obj,  y=y, give.warnings=FALSE)$lambda.est
hold[k,]<-  lambdaTable[1,]
hold2[k,]<-  lambdaTable[6,]
}
matplot( cbind(hold[,2], hold2[,2]),cbind( hold[,4],hold2[,4]),
 xlab="estimated eff. df", ylab="sigma hat", pch=16, col=c("orange3", "green2"), type="p")
yline( sigma, col="grey", lwd=2)

}
}
\keyword{spatial}
