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
\name{QTps}
\alias{QTps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Robust and Quantile smoothing using a thin-plate spline
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function uses the standard thin plate spline function \code{Tps} and a algorithm based on 
psuedo data to compute robust smoothers based on the Huber weight function. By modifying the 
symmetry of the Huber function and changing the scale one can also approximate a quantile 
smoother. This function is experimental in that is not clear how efficient the psuedo-data
algorithm is acheiving convergence to a solution. 
}
\usage{
QTps(x, Y, ..., f.start = NULL, psi.scale = NULL, C = 1, alpha = 0.5, Niterations = 100,
               tolerance = 0.001, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Locations of observations.
}
  \item{Y}{
Observations
}
  \item{\dots}{
Any other arguments to pass to the Tps function.
}
  \item{C}{Scaling for huber robust weighting function. (See below.) Usually it is better to leave this at 1 and 
just modify the scale  \code{psi.scale} according to the size of the residuals. }
  \item{f.start}{
The initial value for the estimated function. If NULL then the constant function at the 
median of \code{Y} will be used. NOTE: This may not be a very good starting vector and a more robust
method would be to use a local robust smoother. 
}
  \item{psi.scale}{
The scale value for the Huber function.  When C=1, this is the point where the Huber weight function will 
change from quadratic to linear. Default is to use the scale \code{.05*mad(Y)}  and \code{C=1} . Very small scales relative to the 
size of the residuals will cause the estimate to approximate a quantile spline. Very large scales will yield the
ordinary least squares spline. 
}
  \item{alpha}{
The quantile that is estimated by the spline. Default is .5 giving a median. Equivalently this parameter controls the slope of the linear wings in the Huber function  \code{2*alpha} for the positive wing and  \code{2*(1-alpha)}
for the negative wing.

}
  \item{Niterations}{
Maximum number of interations of the psuedo data algorithm 
}
  \item{tolerance}{
Convergence criterion based on the relative change in the predicted values of the function estimate. Specifically if the criterion \code{mean(abs(f.hat.new - f.hat))/mean(abs(f.hat))} is less than  \code{tolerance} the iterations re stopped. 
}
  \item{verbose}{
If TRUE intermediate results are printed out.
}
}
\details{
This is an experimentla function that uses the psuedo-value algorithm to compute a class of robust and quantile problems. 

The Thin Plate Spline/ Kriging model through fields is: Y.k= f(x.k) = P(x.k) + Z(x.k) + e.k

 with the goal of estimating the smooth function: f(x)= P(x) + Z(x)

The extension in this function is that e.k can be heavy tailed or have outliers and one would still like a 
robust estimate of f(x). In the quantile approximation (very small scale parameter) f(x) is an estimate of the 
alpha quantile of the conditional distribution of Y given x. 

The algorithm is iterative and involves at each step tapering the residuals in a nonlinear way. 
Let psi.wght  be this tapering function then given an initial estimate of f, f.hat the new data for smoothing is 

\code{ Y.psuedo<- f.hat + psi.scale* psi.wght( Y - f.hat, psi.scale=psi.scale, alpha=alpha)} 
A thin plate spline is now estimated for these data and a new prediction for f is found. This new vector is
used to define new psuedo values. Convergence is achieved when the the subsequent estimates of f.hat do not
change between interations. The advantage of this algorithm is at every step a standard "least squares" thin
plate spline is fit to the psuedo data. Because only the observation vector is changing at each iteration 
Some matrix decompositions need only be found once and the computations at each subsequent iteration are efficient.
At convergence there is some asymptotic theory to suggest that the psuedo data can be fit using the least
squares spline and the standard smoothing techinques are valid. For example one can consider looking at the
cross-validation function for the psuedo-data as a robust version to select a smoothing parameter. This approach
is different from the weighted least squared algorithm used in the \code{qsreg} function. Also \code{qsreg} is only
designed to work with 1-d cubic smoothing splines.   

The  "rho" function indicating the departure from a pure quadratic loss function  has the definition
\preformatted{
qsreg.rho<-function(r, alpha = 0.5, C = 1) 
     temp<- ifelse( r< 0, ((1 - alpha) * r^2)/C ,  (alpha * r^2)/C)
     temp<- ifelse( r >C,   2 * alpha * r - alpha * C, temp)
     temp<- ifelse( r < -C, -2 * (1 - alpha) * r - (1 - alpha) * C, temp)
     temp
} 

The derivative of this function "psi" is 

\preformatted{
 qsreg.psi<- function(r, alpha = 0.5, C = 1)
     temp <- ifelse( r < 0, 2*(1-alpha)* r/C, 2*alpha * r/C )               
     temp <- ifelse( temp > 2*alpha, 2*alpha, temp)
     temp <- ifelse( temp < -2*(1-alpha), -2*(1-alpha), temp)
     temp
}

Note that if C is very small and if alpha = .5 then psi will essentially be 1 for r > 0 and -1 for  r < 0. 
The key feature here is that outside a ceratin range the residual is truncated to a constant value. This is similar
 to the Windsorizing operation in classical robust statistics. 

Another advantage of the psuedo data algotrithm is that at convergence one can just apply all the usual 
generic functions from Tps to the psuedo data fit. For example, predict, surface, print, etc. Some additional
components are added to the Krig/Tps  object, however, for information about the iterations and original data. 
Note that currently these are not reported in the summaries and printing of the output object. 


}
\value{
A \code{Krig} object with additional components:

  \item{yraw}{ Original Y values}
  \item{conv.info}{A vector giving the convergence criterion at each iteration.}
 \item{conv.flag}{If TRUE then convergence criterion was less than the tolerance value.}
 \item{psi.scale}{Scaling factor used for the psi.wght function.}
  \item{value}{Value of alpha.}
}
\references{
Oh, Hee-Seok, Thomas CM Lee, and Douglas W. Nychka. "Fast nonparametric quantile regression with arbitrary smoothing methods." Journal of Computational and Graphical Statistics 20.2 (2011): 510-526.
}
\author{
Doug Nychka
}

\seealso{
qsreg
}
\examples{

data(ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]



# Smoothing fixed at 50 df 
    look1<- QTps( x,y, psi.scale= 15, df= 50)

\dontrun{
# Least squares spline (because scale is so large)
    look2<- QTps( x,y, psi.scale= 100, df= 50)
#
    y.outlier<- y
# add in a huge outlier.
    y.outlier[58]<- 1e5
    look.outlier1<- QTps( x,y.outlier, psi.scale= 15, df= 50)
# least squares spline.
    look.outlier2<- QTps( x,y.outlier, psi.scale=100 , df= 50)
#
    set.panel(2,2)
    surface( look1)
    title("robust spline")
    surface( look2)
    title("least squares spline")
    surface( look.outlier1,  zlim=c(0,250))
    title("robust spline w/outlier") 
    points( rbind(x[58,]), pch="+")
    surface( look.outlier2, zlim=c(0,250))
    title("least squares spline w/outlier")
    points( rbind(x[58,]), pch="+")
    set.panel()
}
# some quantiles
look50 <- QTps( x,y, psi.scale=.5)
look75 <- QTps( x,y,f.start= look50$fitted.values, alpha=.75)


# a simulated example that finds some different quantiles. 
\dontrun{
set.seed(123)
N<- 400
x<- matrix(runif( N), ncol=1)
true.g<- x *(1-x)*2
true.g<- true.g/ mean( abs( true.g))
y<-  true.g + .2*rnorm( N )

look0 <- QTps( x,y, psi.scale=10, df= 15)
look50 <- QTps( x,y, df=15)
look75 <- QTps( x,y,f.start= look50$fitted.values, df=15, alpha=.75)
}

\dontrun{
# this example tests the quantile estimate by Monte Carlo
# by creating many replicate point to increase the sample size. 
# Replicate points are used because the computations for the 
# spline are dominated by the number of unique locations not the 
# total number of points. 
set.seed(123)
N<- 80
M<- 200
x<- matrix( sort(runif( N)), ncol=1)
x<- matrix( rep( x[,1],M), ncol=1)

true.g<- x *(1-x)*2
true.g<- true.g/ mean( abs( true.g))
errors<- .2*(rexp( N*M) -1)
y<- c(matrix(true.g, ncol=M, nrow=N) + .2 *  matrix( errors, ncol=M, nrow=N))

look0 <- QTps( x,y, psi.scale=10, df= 15)
look50 <- QTps( x,y, df=15)
look75 <- QTps( x,y, df=15, alpha=.75)


bplot.xy(x,y, N=25)
xg<- seq(0,1,,200)
lines( xg, predict( look0, x=xg), col="red")
lines( xg, predict( look50, x=xg), col="blue")
lines( xg, predict( look75, x=xg), col="green")
}
\dontrun{
# A comparison with qsreg
  qsreg.fit50<- qsreg(rat.diet$t,rat.diet$con, sc=.5)
  lam<- qsreg.fit50$cv.grid[,1]
  df<- qsreg.fit50$cv.grid[,2]
  M<- length(lam)
  CV<-rep( NA, M)
  M<- length( df)
  fhat.old<- NULL
  for ( k in M:1){
     temp.obj<- QTps(rat.diet$t,rat.diet$con, f.start=fhat.old,  psi.scale=.5, tolerance=1e-6,
     verbose=FALSE, df= df[k])
     cat(k, " ")
     CV[k] <- temp.obj$Qinfo$CV.psuedo
     fhat.old<- temp.obj$fitted.values
  }
  plot( df, CV, type="l", lwd=2)
# psuedo data estimate
  points( qsreg.fit50$cv.grid[,c(5,6)], col="blue")
# alternative CV estimate via reweighted LS
  points( qsreg.fit50$cv.grid[,c(2,3)], col="red")
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
