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

\name{smooth.2d}
\alias{smooth.2d}
\title{
  Kernel smoother for irregular 2-d data
}
\description{
An approximate Nadaraya Watson kernel smoother is obtained by first
discretizing the locations to a grid and then using convolutions to find 
and to apply the kernel weights. The main advantage of this function
is a smoother that avoids explicit looping. 
}
\usage{
smooth.2d(Y, ind = NULL, weight.obj = NULL, setup = FALSE, grid = NULL,
    x = NULL, nrow = 64, ncol = 64, surface = TRUE, cov.function =
gauss.cov, Mwidth = NULL, Nwidth = NULL, ...)
}
\arguments{
\item{Y}{
A vector of 
data to be smoothed   
}
\item{ind}{
Row and column indices that correspond to 
the locations of the data on regular grid. This is most useful when
smoothing the same locations many times. (See also the x argument.)
}
\item{weight.obj }{
An object that 
has the FFT of the convolution kernel and other information ( i.e. the
result from calling this with setup=TRUE).   
}
\item{setup}{
If true creates a list that includes the FFT of the 
convolution kernel. In this case the function will return this
list. Default is false.  
}
\item{grid}{
A list with components x and y 
being equally spaced values that define the grid. Default are integers
1:nrow, 1:ncol. If x is given the ranges will be used to define the grid.  
}
\item{x}{
Actual locations 
of the Y values. Not needed if ind is specified.   
}
\item{nrow}{
 Number of 
points in the horizontal (x) axis of the grid. Not needed if grid is 
specified the default is 64  
}
\item{ncol}{
 Number of points in the vertical (y)
axis of the grid. Not needed if grid list is specified the default is 64 
}
\item{surface}{
If true (the default) a surface object is returned suitable for use by 
image, persp or 
contour functions. If false then just the nrowXncol matrix of smoothed 
values is returned. 
}
\item{cov.function}{
S function describing the kernel function. To be consistent with the other
spatial function this is in the form of a covariance function. The only
assumption is that this be stationary. Default is the (isotropic) Gaussian. 
}
\item{Nwidth}{
The size of the padding  regions of zeroes when computing the 
(exact) convolution of the kernel with the data. The most conservative 
values are 2*nrow and 2*ncol, the default. If the kernel has support of  
say 2L+1 grid points then the padding region need only be of size L+1. 
}
\item{Mwidth}{
See Nwidth.
}
\item{\dots}{
Parameters that are passed to the smoothing kernel. ( e.g. the scale
parameter theta for the
exponential or gaussian)
}
}
\value{
Either a matrix of smoothed values or a surface object.  
The surface object also has a component 'ind' that gives the subscripts of the image 
matrix where the data is present. 
}
\details{
The irregular locations are first discretized to a regular grid ( using
as.image)  
then a 2d- FFT is used to compute a 
Nadaraya-Watson type kernel estimator. Here we take advantage of two
features. The kernel estimator is a convolution and by padding the regular
by zeroes where data is not obsevred one can sum the kernel over irregular
sets of locations. 
A second convolutions to find the normalization of the kernel 
weights.  

The kernel function is specified by an function that should evaluate with
the kernel for two matrices of locations. Assume that the kernel has the
form:  K( u-v) for two locations u and v. The function given as the
argument to cov.function should
have the call myfun( x1,x2) where x1 and x2 are matrices of 2-d locations
if nrow(x1)=m and nrow( x2)=n then this function should return a mXn
matrix where the (i,j) element is K( x1[i,]- x2[j,]).  Optional arguments
that are included in the ... arguments are passed to this function when it
is used. The default kernel is the Gaussian and the argument theta is the
bandwidth. It is easy to write other other kernels, just use 
Exp.cov.simple as
a template. 
}
\examples{
# Normal kernel smooth of the precip data with bandwidth of .5 ( degree) 
#  
look<- smooth.2d( RMprecip$y,  x=RMprecip$x, theta=.25)

# finer resolution used in computing the smooth 
look3<-smooth.2d( RMprecip$y, x=RMprecip$x, theta=.25, nrow=256, 
ncol=256,Nwidth=32,
Mwidth=32) 
# if the width arguments were omitted the padding would create a  
# 512X 512 matrix with the data filled in the upper 256X256 part. 
# with a bandwidth of .25 degrees the normal kernel is essentially zero  
# beyond 32 grid points from its center ( about 6 standard deviations) 
#
# take a look:

#set.panel(2,1)
#image( look3, zlim=c(-8,12))
#points( RMprecip$x, pch=".")  
#image( look, zlim =c(-8,12))
#points( RMprecip$x, pch=".")  


# bandwidth changed to .25, exponential kernel   
look2<- smooth.2d( RMprecip$y, x=RMprecip$x, cov.function=Exp.cov,theta=.25)
# 


 
}
\keyword{smooth}
% docclass is function
% Converted by Sd2Rd version 1.21.
