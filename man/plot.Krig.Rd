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

\name{plot.Krig}
\alias{plot.Krig}
\alias{plot.spatialProcess}
\alias{plot.sreg}

\title{
  Diagnostic and summary plots of a Kriging, spatialProcess or spline object. 
}
\description{
Plots a series of four diagnostic plots that summarize the fit. 
}
\usage{
\method{plot}{Krig}(x, digits=4, which= 1:4,...)
\method{plot}{sreg}(x, digits = 4, which = 1:4, ...)
\method{plot}{spatialProcess}(x, digits=4, which= 1:4,...)


}
\arguments{
\item{x}{ A Krig or an sreg object}

\item{digits}{
Number of significant digits for the RMSE label. 
}
\item{which}{
A vector specifying by number which of the four plots to draw.  
1:4 plots all four. 
}
\item{\dots}{
Optional graphics arguments to pass to each plot. 
}
}
\details{
 This function creates four summary plots of the Krig or sreg object. The
default is to put these on separate pages. However if the screen is
already divided in some other fashion the plots will just be added
according to that scheme. This option is useful to compare to compare
several different model fits. 

The first is a scatterplot of predicted value against observed.

The second plot is "standardized" residuals against predicted value.  
Here we mean that the residuals are divided by the GCV estimate for sigma
and multiplied by the square root of any weights that have been specified.
In the case of a "correlation model" the residuals are also divided by the
marginal standard deviation from this model.

The third plot are the values of the GCV function against the effective
degrees of freedom. When there are replicate points several versions of
the GCV function may be plotted.  GCV function is with respect to the
standardized data if a correlation model is specified. A vertical line
indicates the minimium found.

For \code{Krig} and \code{sreg} objects the fourth plot is a histogram of the standardized residuals. 
For sreg if multiple lambdas are given plotted are boxplots of the  
residuals for each fit. 

For \code{spatialProcess} object the fourth plot is the profile likelihood for the
theta parameter. Points are the actual evaluated log likelihoods and the dashed line is
just a spline interpolation to help with visualization. 

}
\seealso{
Krig, spatialProcess, summary.Krig, Tps, set.panel
}
\examples{

fit<-Krig(ChicagoO3$x, ChicagoO3$y, theta=200)  
# fitting a surface to ozone  
# measurements 
set.panel( 2,2)
plot(fit)

fit<-spatialProcess(ChicagoO3$x, ChicagoO3$y)  
# fitting a spatial process model to ozone  
# measurements 
# Although an example does not make too much sense for only 20 observations!
set.panel( 2,2)
plot(fit)

# fit rat data
fit<-sreg(rat.diet$t,rat.diet$con)
set.panel(2,2)
plot(fit)       

set.panel(1,1) # reset graphics window. 

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
