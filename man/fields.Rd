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

\name{fields}
\alias{fields-package}
\alias{fields}
\title{
fields - tools for spatial data 
}
\description{
 Fields is a collection of programs for curve and function
 fitting with an emphasis on spatial data and spatial statistics. The
 major methods implemented include cubic and thin plate splines, 
 universal
 Kriging and Kriging for large data sets. One main feature is any
 covariance function implemented in R code can be used for spatial prediction. Another important feature is that fields will take advantage of compactly supported covariance functions in a seamless way through 
the spam package. See \code{library( help=fields)} for a listing of all the
fields contents. 

fields stives to have readable and tutorial code. Take a look at the 
source code for \code{Krig} and \code{mKrig} to see how things work 
"under the hood". 
To load fields with the comments retained in the source
use \code{ keep.source = TRUE} in the \code{library} command. 
We also keep the source on-line:
browse the directory 
\url{http://www.image.ucar.edu/~nychka/Fields/Source} for commented source. 
\url{http://www.image.ucar.edu/~nychka/Fields/Help/00Index.html} is a
page for html formatted help files. (If you obtain the source version of the
package (file ends in .gz) the commented source code is the R subdirectory.) 

\strong{Major methods} 
\itemize{ 
	
\item \code{\link{Tps}} Thin Plate spline
regression including GCV and REML estimates for the smoothing parameter. 

\item \code{\link{spatialProcess}} An easy to use method that fits a spatial process model
( e.g. Kriging) but also estimates the key spatial parameters:  nugget variance, sill variance and range by maximum likelihood. Default covariance model is a Matern covariance function.
 
\item \code{\link{Krig}} Spatial process estimation that is the core function of fields.
 
The Krig function allows you to supply a covariance function that is
written in native R code.  See (\code{\link{stationary.cov}}) that includes
several families of covariances and distance metrics including the
Matern and great circle distance. 

\item \code{\link{mKrig}} (micro Krig) are
      \code{\link{fastTps}}
 fast efficient Universal Kriging and spline-like functions, that can take advantage of sparse covariance
functions and thus handle very large numbers of spatial locations. 
 \code{\link{QTps}} A easy to use extension of thin plate splines for quantile and robust surface fitting.

\item \code{\link{mKrig.MLE}} for maximum likelihood estimates of covariance parameters. This function also 
handles replicate fields assumed to be independent realizations at the same locations.  

}

\strong{Other noteworthy functions}
\itemize{

\item \code{\link{vgram}} and \code{\link{vgram.matrix}} find variograms for spatial data (and
with temporal replications.

\item \code{\link{cover.design}} Generates space-filling designs where the distance 
function is expresed in R code.

\item \code{as.image}, \code{image.plot}, \code{drape.plot}, \code{quilt.plot}
\code{add.image}, \code{crop.image}, \code{half.image}, \code{average.image}, 
\code{\link{designer.colors}}, \code{\link{color.scale}}, \code{\link{in.poly}} Many
convenient functions for working with image data and rationally (well,
maybe reasonably) creating and  placing a color scale on an image plot. 
See also  \code{\link{grid.list}} for how fields works with grids and \code{\link{US}}
and \code{\link{world}} for adding a map quickly. 

\item \code{\link{sreg}} \code{\link{splint}}   Fast 1-D smoothing 
splines and interpolating cubic splines.

}


\strong{ Generic functions that support the methods} 

\code{plot} - diagnostic plots of fit \cr
\code{summary}- statistical summary of fit \cr
\code{print}- shorter version of summary \cr
\code{\link{surface}}- graphical display of fitted surface \cr
\code{predict}- evaluation fit at arbitrary points \cr
\code{\link{predictSE}}- prediction standard errors at arbitrary points. \cr
\code{\link{sim.rf}}- Simulate a random fields on a 2-d grid.

\strong{Getting Started}

 Try some of the examples from help files for \code{Tps} or 
\code{spatialProcess}. 

\strong{Graphics tips}

\code{help( fields.hints)}
 gives some R code tricks for setting up common legends and axes. 
And has little to do with this package!

\strong{Testing}
See \code{help(fields.tests)} for testing fields. 

\strong{Some fields datasets}
\itemize{
\item \code{\link{CO2}} Global satelite CO2 concentrations (simulated field)
\item \code{\link{RCMexample}} Regional climate model output
\item \code{\link{lennon}} Image of John Lennon
\item \code{\link{COmonthlyMet}} Monthly mean temperatures and precip for Colorado
\item \code{\link{RMelevation}} Digital elevations for the Rocky Mountain Empire
\item \code{\link{ozone2}}  Daily max 8 hour ozone concentrations for the US midwest 
for summer 1987.
\item \code{\link{PRISMelevation}} Digital elevations for the 
   continental US at approximately 4km resolution
\item \code{\link{NorthAmericanRainfall}} 50+ year average and trend for summer rainfall at 
1700+ stations. 
\item \code{\link{rat.diet}} Small paired study on rat food intake over time.
\item \code{\link{WorldBankCO2}} Demographic and carbon emission data 
 for 75 countries and for 1999. 
}

\strong{DISCLAIMER:}
 The
authors can not guarantee the correctness of any function or program in
this package. 

}

\examples{

# some air quality data, daily surface ozone measurements for the Midwest:
data(ozone2)
x<-ozone2$lon.lat
y<- ozone2$y[16,] # June 18, 1987

# pixel plot of spatial data
quilt.plot( x,y)
US( add=TRUE) # add US map

fit<- Tps(x,y)
# fits a GCV thin plate smoothing spline surface to ozone measurements.
# Hey, it does not get any easier than this!

summary(fit) #diagnostic summary of the fit 
set.panel(2,2)
plot(fit) # four diagnostic plots of fit and residuals.

# quick plot of predicted surface
set.panel()
surface(fit) # contour/image plot of the fitted surface
US( add=TRUE, col="magenta", lwd=2) # US map overlaid
title("Daily max 8 hour ozone in PPB,  June 18th, 1987")


fit2<- spatialProcess( x,y)
# a "Kriging" model. The covariance defaults to a Matern with smoothness 1.0.
# the nugget, sill and range parameters are found by maximum likelihood
# summary, plot, and surface also work for  fit2 !


}
\keyword{datasets}
