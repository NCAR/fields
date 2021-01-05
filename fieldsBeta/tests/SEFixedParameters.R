# fields is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2018
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research,
# PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# tests of SEs from the Omega matrix 
# using Monte Carlo

suppressMessages(library(fields))
options( echo=FALSE)

test.for.zero.flag<- TRUE
X<- ChicagoO3$x
n<- nrow( X)
tau2<- .1
sigma<- 2.0
processCov<-  sigma*Exp.cov( X,X,aRange=50)

cholCov<- chol( processCov + diag( tau2 , n ) )

nreps<- 1e5

set.seed( 111)
# extra random covariate
Z<- matrix( rnorm(n),n)
# the stochastic part 
E<- matrix( rnorm(n*nreps),n,nreps)
# NOTE: all fixed effects set to zero
Y<- t( cholCov)%*%E 
out<- mKrig( X,Y, aRange=50, Z=Z,
             collapseFixedEffect = FALSE,
              lambda=.1/2.0)
testCov<- var( t(out$beta) )

test.for.zero(testCov, sigma*out$Omega, tol=.05  )












