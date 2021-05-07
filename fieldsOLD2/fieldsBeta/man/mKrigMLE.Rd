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
\name{mKrigMLE}
\alias{mKrigMLEJoint}
\alias{mKrigMLEGrid}
\alias{mKrigJointTemp.fn}
\alias{profileCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximizes likelihood for the process marginal variance (sigma) and
 nugget standard deviation (tau) parameters (e.g. lambda) over a
 many covariance models or covariance parameter values.
}
\description{
These function are designed to explore the likelihood surface for
different covariance parameters with the option of maximizing over
tau and sigma. They used the \code{mKrig} base are designed for computational efficiency.
}
\usage{
mKrigMLEGrid(x, y, weights = rep(1, nrow(x)), Z = NULL, 
                 mKrig.args = NULL, cov.function = "stationary.cov",
                 cov.args = NULL, na.rm = TRUE, par.grid = NULL,
                 relative.tolerance = 1e-04,
                 REML = FALSE, GCV=FALSE,
                 optim.args = NULL, 
                 cov.params.start = NULL, verbose = FALSE)

mKrigMLEJoint(x, y, weights = rep(1, nrow(x)), Z = NULL, mKrig.args
                 = NULL, na.rm = TRUE, cov.function = "stationary.cov",
                 cov.args = NULL, cov.params.start = NULL, optim.args =
                 NULL, abstol = 1e-04, parTransform = NULL, 
                 REML = FALSE, GCV=FALSE,
                 verbose = FALSE)

profileCI(obj, parName, CIlevel = 0.95, REML = FALSE) 


mKrigJointTemp.fn(parameters, mKrig.args, cov.args, parTransform,
           parNames, REML = FALSE, GCV=FALSE,  capture.env)
           
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{abstol}{Absolute convergence tolerance used in optim.}
 
 \item{capture.env}{For the ML obective function the frame to save the results of the evaluation. This should be the environment of the function calling optim.}
 
\item{CIlevel}{ Confidence level.}
\item{cov.function}{
The name, a text string, of the covariance function.
}
\item{cov.args}{
The  arguments that would also be included in calls
  to the covariance function to specify the fixed part of 
  the covariance model. This is the form of a list
  E.g.\code{cov.args= list( aRange = 3.5)}
}



\item{cov.params.start}{
A list of initial starts for covariance parameters  to perform likelihood maximization.  The list 
  contains the names of the parameters as well as the values.
  It usually makes sense to optimize over the important lambda  parameter
  (tau^2/ sigma^2)  is most spatial applications and
  so if \code{lambda}  is omitted then the component
  \code{lambda = .5} is added to this list. 
}
 
\item{mKrig.args}{A list of additional parameters to supply to the \code{mKrig} function. E.g.\code{ mKrig.args= list( m=1) }to set the regression function to be a constant function.

\code{mKrig} function that are distinct from the covariance model. 
For example \code{mKrig.args= list( m=1 )} will set the polynomial to be 
just a constant term (degree = m - 1 = 0).
Use  \code{mKrig.args= list( m = 0 )} to omit a fixed model and assume the observations have an expectation of zero. 
}

\item{na.rm}{Remove NAs from data.}

\item{optim.args}{
Additional arguments that would also be included in calls
  to the optim function in joint likelihood maximization.  If 
  \code{NULL}, this will be set to use the "BFGS-" 
  optimization method.  See \code{\link{optim}} for more 
  details.  The default value is: 
  \code{optim.args = list(method = "BFGS", 
             control=list(fnscale = -1, 
                          ndeps = rep(log(1.1), length(cov.params.start)+1), 
                          abstol=1e-04, maxit=20))}
  Note that the first parameter is lambda and the others are 
  the covariance parameters in the order they are given in 
  \code{cov.params.start}.  Also note that the optimization 
  is performed on a transformed scale (based on the function
  \code{parTransform} ), and this should be taken into 
  consideration when passing arguments to \code{optim}.
}
\item{parameters}{The parameter values for evaluate the likelihood.}

 \item{par.grid}{
A list or data frame with components being parameters for
  different covariance models. A typical component is "aRange"
  comprising a vector of scale parameters to try. If par.grid
  is "NULL" then the covariance model is fixed at values that
  are given in \dots.
}

\item{obj}{ List returnerd from \code{mKrigMLEGrid}}
\item{parName}{Name of parameter to find confidence interval.}

\item{parNames}{Names of the parameters to optimize over.}

\item{parTransform}{A function that maps the parameters to a scale
for optimization or
effects the inverse map from the transformed scale into the original values. See below for more details. 
}



\item{relative.tolerance}{
Tolerance used to declare convergence when
  maximizing likelihood over lambda.
}

\item{REML}{If TRUE use REML instead of the full log likelihood.}

\item{GCV}{NOT IMPLEMENTEDm YET! A placeholder for future versions to implement optimization using an approximate cross-validation criterion. }
 


\item{verbose}{ If \code{TRUE} print out interesting intermediate results.
}

\item{weights}{
Precision ( 1/variance) of each observation
}

  \item{x}{

Matrix of unique spatial locations (or in print or surface 
  the returned mKrig object.)
}

  \item{y}{
Vector or matrix of observations at spatial locations, 
  missing values are not allowed! Or in mKrig.coef a new 
  vector of observations. If y is a matrix the columns are 
  assumed to be independent observations vectors generated 
  from the same covariance and measurment error model.
}

\item{Z}{
Linear covariates to be included in fixed part of the 
  model that are distinct from the default low order 
  polynomial in \code{x}
}


}
\details{
The observational model follows the same as that described in the
\code{Krig} function and thus the two primary covariance parameters
for a stationary model are the nugget standard deviation (tau) and
the marginal variance of the process (sigma). It is useful to
reparametrize as \code{sigma} and  \code{ lambda = tau^2/sigma}.
The likelihood can be
maximized analytically over sigma and the parameters in the fixed part
of the model, this estimate of sigma can be substituted back into the
likelihood to give a expression that is just a function of lambda and
the remaining covariance parameters. This operation is called concentrating the likelhood by maximizing over a subset of parameters

For these kind of computations there has to be some device to identify parameters that are fixed and those that are optimized. For \code{mKrigMLEGrid} and \code{mKrigMLEJoint} the list \code{cov.args} should have the fixed parameters.
For example this is how to fix a lambda value in the model. The list
\code{cov.params.start} should be list with all parameters to optimize. The values for each component are use as the starting values. This is how the \link{optim} function works.


For \code{mKrigMLEJoint} the 
 default transformation of the parameters is set up for a log/exp transformation:
\preformatted{
 parTransform <- function(ptemp, inv = FALSE) {
            if (!inv) {
                log(ptemp)
            }
            else {
                exp(ptemp)
            }
        }
}
}

\value{
\strong{\code{mKrigMLEGrid}} returns a list with the components:

\item{summary}{A matrix with each row giving the results of evaluating
the
 likelihood for each covariance model.}

\item{par.grid}{The par.grid argument used. A matrix where rows are the combination of parameters considered.}

\item{call}{The calling arguments to this function.}

\strong{\code{mKrigMLEJoint}} returns a list with components:

\item{summary}{A vector giving the MLEs and the log likelihood at the maximum}

\item{lnLike.eval}{
A table containing information on all likelihood evaluations 
performed in the maximization process.
}
\item{optimResults}{The list returned from the optim function. Note that the parameters may be transformed values. }

\item{par.MLE}{The maximum likelihood estimates.}

\item{parTransform}{The transformation of the parameters used in the optimziation.}

}
\references{
\url{https://github.com/NCAR/Fields}
}
\author{
%%  ~~who you are~~
Douglas W. Nychka, John Paige
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{mKrig}}
\code{\link{Krig}}
\code{\link{stationary.cov}}
\code{\link{optim}}
}

\examples{

\dontrun{
#perform joint likelihood maximization over lambda and aRange. 
# NOTE: optim can get a bad answer with poor initial starts.
  data(ozone2)
  s<- ozone2$lon.lat
  z<- ozone2$y[16,]
  gridList<- list( aRange = seq( .4,1.0,length.out=20),
                 lambda = 10**seq( -1.5,0,length.out=20)
                 )
  par.grid<- make.surface.grid( gridList)
  out<- mKrigMLEGrid( s,z, par.grid=par.grid,
                          cov.args= list(smoothness=1.0,
                                     Covariance="Matern" )
                          )   
  outP<- as.surface( par.grid, out$summary[,"lnProfileLike.FULL"])
  image.plot( outP$x, log10(outP$y),outP$z,
               xlab="aRange", ylab="log10 lambda")
               
}
 
 \dontrun{
  N<- 50
  set.seed(123)
  x<- matrix(runif(2*N), N,2)
  aRange<- .2
  Sigma<-  Matern( rdist(x,x)/aRange , smoothness=1.0)
  Sigma.5<- chol( Sigma)
  tau<- .1
  #  250 independent spatial data sets but a common covariance function 
  #    -- there is little overhead in
  #        MLE across independent realizations and a good test of code validity.
  M<-250
  F.true<- t( Sigma.5) \%*\% matrix( rnorm(N*M), N,M)
  Y<-  F.true +  tau* matrix( rnorm(N*M), N,M)

# find MLE for lambda with grid of ranges 
# and smoothness fixed in Matern                     
 par.grid<- list( aRange= seq( .1,.35,,8))
  obj1b<- mKrigMLEGrid( x,Y,
     cov.args = list(Covariance="Matern", smoothness=1.0), 
        par.grid = par.grid
                    )
  obj1b$summary # take a look
# profile over aRange
  plot( par.grid$aRange, obj1b$summary[,"lnProfileLike.FULL"],
    type="b", log="x")
 }
  \dontrun{
# m=0 is a simple switch to indicate _no_ fixed spatial drift
# (the default and highly recommended  is linear drift, m=2). 
# However, m=0 results in MLEs that are less biased, being the correct model
# -- in fact it nails it !

  obj1a<- mKrigMLEJoint(x,Y, 
                    cov.args=list(Covariance="Matern", smoothness=1.0), 
                    cov.params.start=list(aRange =.5, lambda = .5),
                     mKrig.args= list( m=0))
 
 test.for.zero( obj1a$summary["tau"], tau, tol=.007)
 test.for.zero( obj1a$summary["aRange"], aRange, tol=.015)
 
} 


##########################################################################
# A bootstrap example
# Here is a example of a more efficient (but less robust) bootstrap using 
# mKrigMLEJoint and tuned starting values
##########################################################################
\dontrun{
data( ozone2)
obj<- spatialProcess( ozone2$lon.lat,ozone2$y[16,] )

######### boot strap 
  set.seed(123)
  M<- 250
# create M indepedent copies of the observation vector
  ySynthetic<- simSpatialData( obj, M)
  bootSummary<- NULL
  
 aRangeMLE<- obj$summary["aRange"]
 lambdaMLE<- obj$summary["lambda"]
 
  for(  k in 1:M){
  cat( k, " " )
# here the MLEs are found using the easy top level level wrapper
# see mKrigMLEJoint for a more efficient strategy
  out <- mKrigMLEJoint(obj$x, ySynthetic[,k],
                 weights = obj$weights,
              mKrig.args = obj$mKrig.args,
                 cov.function = obj$cov.function.name,
                cov.args = obj$cov.args, 
        cov.params.start = list( aRange = aRangeMLE,
                                lambda = lambdaMLE)
                      )
  newSummary<- out$summary
  bootSummary<- rbind( bootSummary, newSummary)
  }
  
  cat(  " ", fill=TRUE )
  
  obj$summary
  stats( bootSummary)
  
}
\dontrun{
#perform joint likelihood maximization over lambda, aRange, and smoothness.  
#note: finding smoothness is not a robust optimiztion 
#      can get a bad answer with poor initial guesses.
obj2<- mKrigMLEJoint(x,Y, 
                      cov.args=list(Covariance="Matern"), 
                      cov.params.start=list( aRange = .18,
                                         smoothness = 1.1,
                                             lambda = .08),
                       )

#look at lnLikelihood  evaluations
obj2$summary
#compare to REML
obj3<- mKrigMLEJoint(x,Y, 
                      cov.args=list(Covariance="Matern"), 
                      cov.params.start=list(aRange = .18, 
                                       smoothness = 1.1,
                                           lambda = .08),
                       , REML=TRUE)
obj3$summary                      
}
\dontrun{
#look at lnLikelihood  evaluations

# check convergence of MLE to true fit with no fixed part
# 
obj4<- mKrigMLEJoint(x,Y, 
                      mKrig.args= list( m=0),
                      cov.args=list(Covariance="Matern", smoothness=1),
                      cov.params.start=list(aRange=.2, lambda=.1),
                       REML=TRUE)
#look at lnLikelihood  evaluations
obj4$summary
# nails it!
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{spatial}
