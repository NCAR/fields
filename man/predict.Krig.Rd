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

\name{predict.Krig}
\alias{predict.Krig}
\alias{predict.Tps}
\alias{predictDerivative.Krig}
\alias{predict.fastTps}
\title{
  Evaluation of Krig spatial process estimate.  
}
\description{
Provides predictions from the Krig spatial process estimate at arbitrary
points, new data (Y) or other values of the smoothing parameter (lambda)
including a GCV estimate. 
}
\usage{
\method{predict}{Krig}(
object, x = NULL, Z = NULL, drop.Z = FALSE, just.fixed
                 = FALSE, lambda = NA, df = NA, model = NA,
                 eval.correlation.model = TRUE, y = NULL, yM = NULL,
                 verbose = FALSE, ...)
predictDerivative.Krig(object, x = NULL,  verbose = FALSE,...)

\method{predict}{Tps}(object, ... )

\method{predict}{fastTps}(object, xnew = NULL, grid.list = NULL, ynew = NULL,
                 derivative = 0, Z = NULL, drop.Z = FALSE, just.fixed =
                 FALSE, xy = c(1, 2), ...)

}
\arguments{

\item{derivative}{The degree of the derivative to be evauated. Default is
0 (evaluate the function itself), 1 is supported by some covariance functions,
Higher derivatives are not supported in this version and for mKrig.}

\item{df}{
Effective degrees of freedom for the predicted surface. This can be used
in place of lambda ( see the function Krig.df.to.lambda)
}

\item{eval.correlation.model}{
If true ( the default) will multiply the predicted function by marginal
sd's
and add the mean function. This usually what one wants. If false will
return predicted surface in the standardized scale. The main use of this
option is a call from Krig to find MLE's of rho and sigma2     
}

\item{grid.list}{A \code{grid.list} specfiying a grid of locations to
evaluate the fitted surface.}


\item{just.fixed}{ Only fixed part of model is evaluated}


\item{lambda}{ Smoothing parameter. If omitted, out\$lambda will be
used.  (See also df and gcv arguments) } \item{model}{ Generic
argument that may be used to pass a different lambda.  }
\item{object}{ Fit object from the Krig, Tps, mKrig, or fastTps
functions.  }

 \item{verbose}{ Print out all kinds of intermediate stuff for
debugging }

\item{xy}{The column positions that locate the x and y variables for
evaluating on a grid.  This is mainly useful if the surface has more
than 2 dimensions.}

\item{y}{ Evaluate the estimate using the new data vector y (in the
same order as the old data). This is equivalent to recomputing the
Krig object with this new data but is more efficient because many
pieces can be reused. Note that the x values are assumed to be the
same.  } \item{x}{ Matrix of x values on which to evaluate the kriging
surface.  If omitted, the data x values, i.e. out\$x will be used.  }

\item{xnew}{Same as x above.}

\item{ynew}{Same as y above.}  

\item{yM}{ If not NULL evaluate the
estimate using this vector as the replicate mean data. That is, assume
the full data has been collapsed into replicate means in the same
order as xM. The replicate weights are assumed to be the same as the
original data. (weightsM) }


\item{Z}{ Vector/Matrix of additional covariates to be included in
fixed part of spatial model} \item{drop.Z}{ If TRUE only spatial fixed
part of model is evaluated.  i.e. Z covariates are not used.  }

\item{\dots}{Other arguments passed to covariance function. In the case of
\code{fastTps} these are the same arguments as \code{predict.mKrig}.
This argument is usually not needed. 
}

}
\value{
Vector of predicted responses or a matrix of the partial derivatives. 
}
\details{
 The main goal in this function is to reuse the Krig object to rapidly 
evaluate different estimates. Thus there is flexibility in changing the 
value of lambda and also the independent data without having to 
recompute the matrices associated with the Krig object. The reason this 
is possible is that most on the calculations depend on the observed 
locations not on lambda or the observed data. Note the version for 
evaluating partial derivatives does not provide the same flexibility as 
\code{predict.Krig} and makes some assumptions about the null model 
(as a low order polynomial) and can not handle the correlation model form.

}
\seealso{
Krig, predictSurface gcv.Krig 
}

\examples{
  Krig(ChicagoO3$x,ChicagoO3$y, theta=50) ->fit
  predict( fit) # gives predicted values at data points should agree with fitted.values
                #  in fit object 

# predict at the coordinate (-5,10)
  x0<- cbind( -5,10) # has to be a  1X2 matrix
  predict( fit,x= x0)

# redoing predictions at data locations:
   predict( fit, x=ChicagoO3$x)

# only the fixed part of the model
  predict( fit, just.fixed=TRUE) 

# evaluating estimate at a grid of points 
  grid<- make.surface.grid( list( seq( -40,40,,15), seq( -40,40,,15)))
  look<- predict(fit,grid) # evaluate on a grid of points

# some useful graphing functions for these gridded predicted values
  out.p<- as.surface( grid, look) # reformat into $x $y $z image-type object
  contour( out.p) 

# see also the functions predictSurface and surface 
# for functions that combine these steps 
   

# refit with 10 degrees of freedom in surface
  look<- predict(fit,grid, df=15)
# refit with random data 
  look<- predict( fit, grid, y= rnorm( 20))


# finding partial derivatives of the estimate
#
# find the partial derivatives at observation locations
# returned object is a two column matrix. 
# this does not make sense for the exponential covariance
# but can illustrate this with a thin plate spline with
# a high enough order ( i.e. need m=3 or greater)
# 
  data(ozone2)
# the 16th day of this ozone spatial dataset
  fit0<- Tps( ozone2$lon.lat, ozone2$y[16,], m=3)
  look1<- predictDerivative.Krig( fit0)
# for extra credit compare this to
  look2<- predictDerivative.Krig( fit0, x=ozone2$lon.lat)  
# (why are there more values in look2) 


}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
