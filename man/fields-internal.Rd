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

\name{fields internal }
\alias{[.spatialDesign}
\alias{COR}
\alias{D4transform.image}
\alias{Krig.df.to.lambda}
\alias{Krig.fdf}
\alias{Krig.fgcv}
\alias{Krig.fgcv.model}
\alias{Krig.fgcv.one}

\alias{Krig.flplike}
\alias{Krig.fs2hat}
\alias{Krig.ftrace}
\alias{Krig.parameters}
\alias{Krig.updateY}
\alias{Krig.which.lambda}
\alias{Krig.ynew}
\alias{bisection.search}
\alias{cat.matrix}
\alias{cat.to.list}
\alias{ceiling2}
\alias{describe}
\alias{dyadic.2check}
\alias{dyadic.check}
\alias{double.exp}
\alias{Exp.earth.cov}
\alias{fast.1way}
\alias{find.upcross}
\alias{gauss.cov}
\alias{golden.section.search}
\alias{imageplot.info}
\alias{imagePlotInfo}
\alias{imageplot.setup}
\alias{krig.image.parameters}

\alias{makeSimulationGrid}
\alias{makeSimulationGrid2}
\alias{makePredictionPoints}
\alias{multWendlandGrid}

\alias{minimax.crit}
\alias{plot.krig.image}
\alias{plot.sim.krig.image}
\alias{plot.spatialDesign}
\alias{predict.interp.surface}
\alias{predict.krig.image}

\alias{predict.surface.default}

\alias{print.krig.image}
\alias{print.spatialDesign}
\alias{print.sreg}
\alias{print.summary.Krig}
\alias{print.summary.spatialProcess}
\alias{print.summary.krig.image}
\alias{print.summarySpatialDesign}
\alias{print.summary.sreg}
\alias{printGCVWarnings}
\alias{qr.q2ty}
\alias{qr.yq2}
\alias{plot.qsreg}
\alias{predict.qsreg}
\alias{print.qsreg}
\alias{qsreg.fit}
\alias{qsreg.psi}
\alias{qsreg.rho}
\alias{qsreg.psi.OLD}
\alias{qsreg.rho.OLD}
\alias{qsreg.trace}
\alias{quickPrint}
\alias{summary.qsreg}
\alias{radbas.constant}
\alias{sim.krig.image}
\alias{sreg.df.to.lambda}
\alias{sreg.fdf}
\alias{sreg.fgcv}
\alias{sreg.fgcv.model}
\alias{sreg.fgcv.one}
\alias{sreg.fit}
\alias{sreg.fs2hat}
\alias{sreg.trace}
\alias{stats.sim.krig.image}
\alias{summaryGCV.Krig}
\alias{summaryGCV.sreg}
\alias{summary.krig.image}
\alias{summary.spatialDesign}
\alias{summary.sreg}
\alias{surface}
\alias{surface.default}
\alias{surface.krig.image}
\alias{unscale}
\alias{world.dat}
\alias{compactTOMatOLD}

\alias{MLESpatialProcess.fast}

\title{
  Fields internal and secondary functions 
}
\description{
Listed below are supporting functions for the major methods in fields. 
}
\usage{
\method{[}{spatialDesign}(x, ...)

Krig.df.to.lambda(df, D, guess = 1, tol = 1e-05) 
Krig.fdf (llam, info) 
Krig.fgcv (lam, obj) 
Krig.fgcv.model (lam, obj) 
Krig.fgcv.one (lam, obj) 
Krig.flplike (lambda, obj)
Krig.fs2hat (lam, obj) 
Krig.ftrace (lam, D) 
Krig.parameters (obj, mle.calc=obj$mle.calc)
Krig.updateY (out, Y, verbose = FALSE, yM=NA) 
Krig.which.lambda(out)
Krig.ynew (out, y=NULL, yM=NULL ) 
 



bisection.search (x1, x2, f, tol = 1e-07, niter = 25, f.extra = 
        NA, upcross.level = 0) 
 
cat.matrix (mat, digits = 8) 
 
cat.to.list (x, a) 
 
ceiling2 (m) 
 
describe (x) 

double.exp(x)

dyadic.2check( m,n,cut.p=2)
dyadic.check( n,cut.p=2) 

Exp.earth.cov (x1, x2, theta = 1) 
 
fast.1way  (lev, y, w = rep(1, length(y))) 
 
find.upcross (fun, fun.info, upcross.level = 0, guess = 1, tol = 
1e-05) 
 
gauss.cov (...) 


golden.section.search (ax, bx, cx, f, niter = 25, f.extra = NA, 
     tol = 1e-05, gridx=NA) 
 
imagePlotInfo (...,breaks, nlevel) 
imageplot.info(...)
imageplot.setup(x, add=FALSE,  legend.shrink = 0.9, legend.width = 1, 
horizontal = FALSE, legend.mar=NULL, bigplot = NULL, smallplot = NULL,...)


makeSimulationGrid(mKrigObject, predictionPoints, nx, ny, nxSimulation,
                 nySimulation, gridRefinement, gridExpansion)
makeSimulationGrid2 (fastTpsObject, predictionPointsList,
                               gridRefinement, gridExpansion)


  
minimax.crit (obj, des = TRUE, R) 
 

\method{plot}{spatialDesign}(x,...) 
 
\method{predict}{interp.surface}(object, loc,...)  

\method{predict}{surface}(object, ...) 
\method{predict}{surface.default}(object, ...) 
 

\method{print}{spatialDesign} (x,...) 
\method{print}{sreg}(x, ...) 
\method{print}{summary.Krig} (x, ...) 
\method{print}{summary.spatialProcess} (x, ...)
\method{print}{summarySpatialDesign} (x, digits = 4,...) 
\method{print}{summary.sreg} (x, ...) 
 
printGCVWarnings( Table,  method = "all")
 
makePredictionPoints(mKrigObject, nx, ny, predictionPointsList)

multWendlandGrid( grid.list,center, delta, coef, xy = c(1, 2))

qr.q2ty (qr, y) 
 
qr.yq2 (qr, y) 
\method{plot}{qsreg}(x, pch = "*", main = NA,...) 
\method{predict}{qsreg}(object, x, derivative = 0, model = object$ind.cv.ps,...)
\method{print}{qsreg} (x, ...) 
qsreg.fit (x, y, lam, maxit = 50, maxit.cv = 10, 
           tol = 1e-04, offset = 0, sc = sqrt(var(y)) * 1e-07, alpha = 0.5, 
           wt = rep(1, length(x)), cost = 1)
qsreg.psi( r,alpha=.5,C=1) 
qsreg.rho( r,alpha=.5,C=1)
qsreg.trace(x, y, lam, maxit = 50, maxit.cv = 10, tol = 1e-04, 
             offset = 0, sc = sqrt(var(y)) * 1e-07, alpha = 0.5, 
             wt = rep(1, length(x)), cost = 1)
qsreg.rho.OLD(r, alpha = 0.5, C = 1)
qsreg.psi.OLD(r, alpha = 0.5, C = 1)

quickPrint(obj, max.values = 10)



radbas.constant (m, d) 
 
sreg.df.to.lambda (df, x, wt, guess = 1, tol = 1e-05) 
sreg.fdf (h, info) 
sreg.fgcv (lam, obj) 
sreg.fgcv.model (lam, obj) 
sreg.fgcv.one (lam, obj) 
sreg.fit (lam, obj, verbose=FALSE) 
sreg.fs2hat (lam, obj) 
sreg.trace (h, info) 

summaryGCV.Krig(object, lambda, cost = 1, verbose = FALSE, 
                offset = 0, y = NULL, ...)
summaryGCV.sreg (object, lambda, cost = 1, nstep.cv = 20, 
        offset = 0, verbose = TRUE,...) 

\method{summary}{qsreg} (object, ...)     
\method{summary}{spatialDesign} (object, digits = 4, ...) 
\method{summary}{sreg} (object, digits = 4, ...) 

surface(object , ...) 
\method{surface}{default} (object, ...) 


unscale (x, x.center, x.scale) 

MLESpatialProcess.fast(x, y, lambda.start=.5, theta.start = NULL, 
  cov.function = "stationary.cov", 
  cov.args = list(Covariance = "Matern", smoothness = 1), 
  Distance = "rdist", verbose=FALSE, optim.args=NULL, ...)

}
\keyword{internal}
% docclass is function
% Converted by Sd2Rd version 1.21.
