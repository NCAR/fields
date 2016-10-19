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
\name{mKrigMLE}
\alias{mKrigMLEJoint}
\alias{mKrigMLEGrid}
\alias{fastTpsMLE}
\alias{mKrigJointTemp.fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximizes likelihood for the process marginal variance (rho) and
 nugget standard deviation (sigma) parameters (e.g. lambda) over a
 many covariance models or covariance parameter values.
}
\description{
These function are designed to explore the likelihood surface for
different covariance parameters with the option of maximizing over
sigma and rho. They used the \code{mKrig} base are designed for computational efficiency.
}
\usage{
mKrigMLEGrid(x, y, weights = rep(1, nrow(x)), Z = NULL, mKrig.args = NULL, 
    cov.fun = "stationary.cov", cov.args = NULL, na.rm = TRUE, 
    par.grid = NULL, lambda = NULL, lambda.profile = TRUE,
    relative.tolerance = 1e-04, 
    REML = FALSE, verbose = FALSE)

mKrigMLEJoint(x, y, weights = rep(1, nrow(x)), Z = NULL, mKrig.args
                 = NULL, na.rm = TRUE, cov.fun = "stationary.cov",
                 cov.args = NULL, lambda.start = 0.5, cov.params.start
                 = NULL, optim.args = NULL, abstol = 1e-04,
                 parTransform = NULL, REML = FALSE, verbose = FALSE)


fastTpsMLE(x, y, weights = rep(1, nrow(x)), Z = NULL, ...,
                 par.grid=NULL, theta, lambda = NULL, lambda.profile = TRUE,
                 verbose = FALSE, relative.tolerance = 1e-04)

mKrigJointTemp.fn(parameters, mKrig.args, cov.args, parTransform,
           parNames, REML = FALSE, capture.env)
           
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
Matrix of unique spatial locations (or in print or surface 
  the returned mKrig object.)
}
  \item{y}{
%%     ~~Describe \code{y} here~~
Vector or matrix of observations at spatial locations, 
  missing values are not allowed! Or in mKrig.coef a new 
  vector of observations. If y is a matrix the columns are 
  assumed to be independent observations vectors generated 
  from the same covariance and measurment error model.
}
\item{mKrig.args}{A list of additional parameters to supply to the base 
\code{mKrig} function that are distinct from the covariance model. 
For example \code{mKrig.args= list( m=1 )} will set the polynomial to be 
just a constant term (degree = m -1 = 0).
}
\item{cov.fun}{
The name, a text string, of the covariance function.
}
\item{cov.args}{
Additional arguments that would also be included in calls
  to the covariance function to specify the fixed part of 
  the covariance model.
}
  \item{weights}{
%%     ~~Describe \code{weights} here~~
Precision ( 1/variance) of each observation
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
Linear covariates to be included in fixed part of the 
  model that are distinct from the default low order 
  polynomial in \code{x}
}
  \item{par.grid}{
A list or data frame with components being parameters for
  different covariance models. A typical component is "theta"
  comprising a vector of scale parameters to try. If par.grid
  is "NULL" then the covariance model is fixed at values that
  are given in \dots.
}
\item{cov.params.start}{
A list of initial starts for covariance parameters over which 
  the user wishes to perform likelihood maximization.  The list 
  contains the names of the parameters as well as the values.
}
  \item{lambda}{
If \code{lambda.profile=FALSE} the values of lambda to evaluate
  the likelihood if "TRUE" the starting values for the
  optimization. If lambda is NA then the optimum value from
  previous search is used as the starting value. If lambda is
  NA and it is the first value the starting value defaults to
  1.0.
}
  \item{lambda.start}{
The initial guess for lambda in the joint log-likelihood 
  maximization process
}
  \item{lambda.profile}{
  If \code{TRUE} maximize likelihood over lambda.
}
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
  is performed on a transformed scale (based on the function \code{parTransform} ), and this should be taken into 
  consideration when passing arguments to \code{optim}.
}
\item{parTransform}{A function that maps the parameters to a scale for optimization or
effects the inverse map from the transformed scale into the original values. See below for more details. 
}

  \item{REML}{Currently using  REML is not implemented.}
  
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
If \code{TRUE} print out interesting intermediate results.
}
  \item{relative.tolerance}{
Tolerance used to declare convergence when
  maximizing likelihood over lambda.
}
  \item{abstol}{Absolute convergence tolerance used in optim.}
  
  \item{theta}{Range parameter for compact Wendland covariance. (see
  fastTps)}
}
\details{
The observational model follows the same as that described in the
\code{Krig} function and thus the two primary covariance parameters
for a stationary model are the nugget standard deviation (sigma) and
the marginal variance of the process (rho). It is useful to
reparametrize as rho and\ lambda= sigma^2/rho. The likelihood can be
maximized analytically over rho and the parameters in the fixed part
of the model the estimate of rho can be substituted back into the
likelihood to give a expression that is just a function of lambda and
the remaining covariance parameters. It is this expression that is
then maximized numerically over lambda when \code{ lambda.profile =
TRUE}.

Note that \code{fastTpsMLE} is a convenient variant of this more general
version to use directly with fastTps, and \code{mKrigMLEJoint} is 
similar to \code{mKrigMLEGrid}, except it uses the \code{optim} function 
to optimize over the specified covariance parameters and lambda jointly 
rather than optimizing on a grid.  Unlike \code{mKrigMLEJoint}, it returns 
an mKrig object.

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

\item{summary}{A matrix giving the results for evaluating the
 likelihood for each covariance model.}

\item{par.grid}{The par.grid argument used.}

\item{cov.args.MLE}{The list of covariance arguments (except for
lambda) that have the largest likelihood over the list covariance
models. NOTE: To fit the surface at the largest likelihood among those tried
\code{ do.call( "mKrig", c(obj$mKrig.args,
obj$cov.args.MLE,list(lambda=obj$lambda.opt)) )} where \code{obj} is
the list returned by this function.}

\item{call}{The calling arguments to this function.}

\strong{\code{mKrigMLEJoint}} returns a list with components:

\item{summary}{A vector giving the MLEs and the log likelihood at the maximum}

\item{lnLike.eval}{
A table containing information on all likelihood evaluations 
performed in the maximization process.
}
\item{optimResults}{The list returned from the optim function.}

\item{par.MLE}{The maximum likelihood estimates.}

\item{parTransform}{The transformation of the parameters used in the optimziation.}

}
\references{
%% ~put references to the literature/web site here ~
http://cran.r-project.org/web/packages/fields/fields.pdf
http://www.image.ucar.edu/~nychka/Fields/
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
#perform joint likelihood maximization over lambda and theta. 
#optim can get a bad answer with poor initial starts.
data(ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]
obj<- mKrigMLEJoint(x,y, 
                      cov.args=list(Covariance="Matern", smoothness=1.0), 
                      cov.params.start=list(theta=.2), lambda.start=.1)

# check  lnLikeihood evaluations that were culled from optim
# these are in obj$lnLike.eval
# funny ranges are set to avoid outlining and very low likelihood value

quilt.plot( log10(cbind(obj$lnLike.eval[,1:2])), obj$lnLike.eval[,5],
             xlim=c(-1.2,-.40), ylim=c( -1,1), zlim=c( -625, -610))

# some synthetic data
  N<- 100
  set.seed(123)
  x<- matrix(runif(2*N), N,2)
  theta<- .2
  Sigma<-  Matern( rdist(x,x)/theta , smoothness=1.0)
  Sigma.5<- chol( Sigma)
  sigma<- .1
  M<-50 #  fifty independent spatial data sets -- little overhead in MLE across
  # indedpendent realizations.
  #F.true<- t( Sigma.5)%*% matrix( rnorm(N*M), N,M)  
  F.true<- t( Sigma.5)\%*\% matrix( rnorm(N*M), N,M)
  Y<-  F.true +  sigma* matrix( rnorm(N*M), N,M)
# find MLE for lambda with grid of ranges 
# and smoothness fixed in Matern

 par.grid<- list( theta= seq( .05,.3,,20))
  obj<- mKrigMLEGrid( x,Y, cov.args= list(Covariance="Matern", smoothness=1.0), 
        par.grid = par.grid
                    )
  obj$summary # take a look
# profile over theta
  plot( par.grid$theta, obj$summary[,"lnProfileLike.FULL"], type="b")

\dontrun{
#perform joint likelihood maximization over lambda, theta, and smoothness.  
#optim can get a bad answer with poor initial guesses.
obj2<- mKrigMLEJoint(x,Y, 
                      cov.args=list(Covariance="Matern", smoothness=1), 
                      cov.params.start=list(theta=.2),
                       lambda.start=.1)

#look at lnLikelihood  evaluations
obj2$summary
#compare to REML
obj3<- mKrigMLEJoint(x,Y, 
                      cov.args=list(Covariance="Matern", smoothness=1), 
                      cov.params.start=list(theta=.2),
                       lambda.start=.1, REML=TRUE)
}
\dontrun{
#look at lnLikelihood  evaluations
obj3$summary
# check convergence of MLE to true fit with no fixed part
# 
obj4<- mKrigMLEJoint(x,Y, 
                      mKrig.args= list( m=0),
                      cov.args=list(Covariance="Matern", smoothness=1), 
                      cov.params.start=list(theta=.2),
                       lambda.start=.1, REML=TRUE)
#look at lnLikelihood  evaluations
obj4$summary
# nails it!
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{spatial}
