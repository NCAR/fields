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

\name{Krig.Amatrix}
\alias{Krig.Amatrix}
\title{
  Smoother (or "hat") matrix relating predicted values to the dependent 
(Y) values. 
}
\description{
For a fixed value of the smoothing parameter or the covariance 
function some nonparametric curve estimates are linear functions of 
the observed data. This is a intermediate  level function that 
computes the linear weights to be applied to the observations to 
estimate the curve at a particular point. For example the predicted 
values can be represented as Ay where A is an N X N matrix of 
coefficients and Y is the vector of observed dependent variables. For 
linear smoothers the matrix A may depend on the smoothing parameter ( 
or covariance function and the independent variables (X) but NOT on Y.    
}
\usage{
Krig.Amatrix(object, x0 = object$x, lambda=NULL,
      eval.correlation.model = FALSE,...)
}
\arguments{
Output object from fitting a data set using a FIELD regression
method. 
Currently this is supported only for Krig ( and Tps) functions. 

\item{object}{
A Krig object produced by the Krig ( or Tps) function. 
}
\item{x0}{
Locations for prediction default is the observation locations. 
}
\item{lambda}{
Value of the smoothing parameter.
}
\item{eval.correlation.model}{This applies to a correlation model 
where the observations have been standardized -- e.g. 
  y standardized  = (yraw - mean) / (standard deviation). If TRUE the 
prediction in the correlation scale is transformed by the standard
deviation and mean to give a  prediction in the raw scale. 
If FALSE predictions are left in the correlation  scale.}

\item{\dots}{ Other arguments that can used by predict.Krig.}

}
\value{
A matrix where the number of rows is equal to the number of predicted points 
and the number of columns is equal to the length of the Y vector.  
}
\details{
The main use of this function is in finding prediction standard errors. 

For the Krig ( and Tps) functions the A matrix is constructed based on the
representation of the estimate as a generalized ridge regression.  The
matrix expressions are explained in the references from the FIELDS manual.
For linear regression the matrix that gives predicted values is often
referred to as the "hat" matrix and is useful for regression diagnostics.  
For smoothing problems the effective number of parameters in the fit is
usually taken to be the trace of the A matrix. Note that while the A
matrix is usually constructed to predict the estimated curve at the data
points Amatrix.Krig does not have such restrictions. This 
is possible
because any value of the estimated curve will be a linear function of Y.  

The actual calculation in this function is simple. It invovles 
loop through the unit vectors at each observation and computation of the
prediction for each of these delta functions. This approach makes it easy to 
handle different options such as including covariates. 

}
\section{References}{
Nychka (2000) "Spatial process estimates as smoothers." 
}
\seealso{
Krig, Tps, predict.Krig 
}
\examples{
# Compute the A matrix or "hat" matrix for a thin plate spline 
# check that this gives the same predicted values  
tps.out<-Tps( ChicagoO3$x, ChicagoO3$y)
A<-Krig.Amatrix( tps.out, ChicagoO3$x)
test<- A\%*\%ChicagoO3$y 
# now compare this to predict( tps.out) or tps.out$fitted.values 
#                    they should be the same 
stats( test- tps.out$fitted.values)
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
