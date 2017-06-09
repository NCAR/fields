%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2017
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

\name{image.smooth}
\alias{image.smooth}
\alias{setup.image.smooth}
\title{
  Kernel smoother for irregular 2-d data 
}
\description{
Takes an image matrix and applies a kernel smoother to it. Missing values
are handled using the Nadaraya/Watson normalization of the kernel. 
}
\usage{
\method{image}{smooth}(x, wght = NULL, dx = 1, dy = 1,
    kernel.function = double.exp,
    theta = 1, grid = NULL, tol = 1e-08, xwidth = NULL, ywidth = NULL,
    weights = NULL,...)

setup.image.smooth(nrow = 64, ncol = 64, dx = 1, dy = 1,
                   kernel.function = double.exp,
                   theta = 1, xwidth = nrow * dx, ywidth = ncol * dx, lambda=NULL, ...)}
\arguments{
\item{x}{
A matrix image. Missing values can be indicated by NAs. }
\item{wght}{
FFT of smoothing kernel. If this is NULL the default is to compute this
object. }
\item{grid}{
A list with x and y components. Each are equally spaced and define the rectangular. ( see grid.list)}
\item{dx}{
Grid spacing in x direction }
\item{dy}{
Grid spacing in x direction }
\item{kernel.function}{
 An R function that takes as its argument the \emph{squared} distance 
between two points divided by the bandwidth. The default is exp(
-abs(x)) yielding a normal kernel}
\item{theta}{the bandwidth or scale parameter.}
\item{xwidth}{
 Amount of zero padding in horizontal dimension in units of the grid spacing.
If NULL the default value is equal to the width of the image the most 
conservative value but possibly inefficient for computation.
Set this equal to zero to get periodic wrapping of the smoother. This is
useful to smooth a Mercator map projection. }
\item{ywidth}{
 Same as xwidth but for the vertical dimension. }
\item{weights}{
Weights to apply when smoothing.}
\item{tol}{
Tolerance for the weights of the N-W kernel. This avoids kernel
estimates that are "far" away from data. Grid points with weights
less than tol are set to NA.}
\item{nrow}{X dimension of image in setting up smoother weights}
\item{ncol}{Y dimension of image}
\item{lambda}{Smoothing parameter if smoother is interpreted in a spline-like
way.}
\item{\dots}{
 Other arguments to be passed to the kernel function}
}

\value{
The smoothed image in R image format. ( A list with components x, y
and z.)  \code{setup.image.smooth} returns a list with components W a
matrix being the FFT of the kernel, dx, dy, xwidth and ywidth.}

\details{
 The function works by taking convolutions using an FFT. The missing
pixels are taken into account and the kernel smoothing is correctly
normalized for the edge effects following the classical Nadaraya-Watson
estimator. For this reason the kernel doe snot have to be a desity as it
is automatically normalized when the kernel weight function is found for
the data.  If the kernel has limited support then the width arguments
can be set to reduce the amount of computation. (See example below.) 
For multiple smoothing compute the fft of the kernel just once using
\code{setup.image.smooth} and pass this as the wght argument to
image.smooth.  this will save an FFT in computations. 
}
\seealso{ as.image, sim.rf, image.plot}

\examples{
# first convert precip data to the 128X128 discretized image format ( with 
# missing  values to indicate where data is not observed) 
# 
out<- as.image( RMprecip$y, x= RMprecip$x, nx=128, ny=128) 
# out$z is the image matrix 

dx<- out$x[2]- out$x[1] 
dy<-  out$y[2] - out$y[1] 

#  
# grid scale in degrees and choose kernel bandwidth to be .25 degrees. 

look<- image.smooth( out, theta= .25) 

image.plot(look) 
points( RMprecip$x)
US( add=TRUE, col="grey", lwd=2)

# to save on computation, decrease the padding with zeroes 
# only pad 32 grid points around the margins ofthe image. 

look<- image.smooth(out$z, dx=dx, dy=dy, theta= .25, xwidth=32*dx,ywidth=32*dy) 

# the range of these data is ~ 10 degrees  and so 
# with a padding of 32 grid points  32*( 10/128) =  2.5 
# about 10 standard deviations of the normal kernel so there is still 
# lots of room for padding  
# a minimal choice might be  xwidth = 4*(.25)= 1  4 SD for the normal kernel
# creating weighting object outside the call  
# this is useful when one wants to smooth different data sets but on the 
# same grid with the same kernel function 
# 

#
#  random fields from smoothing white noise with this filter.
#
set.seed(123)
test.image<- matrix( rnorm(128**2),128,128)
dx<- .1
dy<- .8

wght<- setup.image.smooth( nrow=128, ncol=128,  dx=dx, dy=dy,
             theta=.25, xwidth=2.5, ywidth=2.5)
#
look<- image.smooth( test.image, dx=dx, dy=dy, wght)

# NOTE:   this is the same as using 
#
#     image.smooth( test.image , 128,128), xwidth=2.5,
#                        ywidth=2.5, dx=dx,dy=dy, theta=.25)
#
#   but the call to image.smooth is faster because fft of kernel
#   has been precomputed.



# periodic smoothing in the horizontal dimension

look<- image.smooth( test.image , xwidth=1.5,
                        ywidth=2.5, dx=dx,dy=dy, theta=1.5)
look2<- image.smooth( test.image , xwidth=0,
                        ywidth=2.5, dx=dx,dy=dy, theta=1.5)
# compare these two
set.panel( 1,2)
image.plot( look, legend.mar=7.1)
title("free boundaries")
image.plot( look2, legend.mar=7.1) # look for periodic continuity at edges!
title("periodic boundary in horizontal")
set.panel(1,1)

}
\keyword{smooth}
% docclass is function
% Converted by Sd2Rd version 1.21.
