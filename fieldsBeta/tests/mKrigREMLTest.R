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


suppressMessages(library(fields))
options( echo=FALSE)
set.seed( 123)
x<- matrix( runif( 20),10, 2)
y<- rnorm(10)


lambda<- .1
aRange<- .2
out<- mKrig( x,y, aRange= aRange, lambda=lambda)

test.for.zero( out$lnDetOmega,
2*log( prod(diag(chol(out$Omega))))
)
    
Mc<-   exp( -rdist( x,x)/aRange) + lambda* diag( 1,10)
OmegaTest<- solve(t(out$Tmatrix)%*%solve( Mc)%*% out$Tmatrix)

test.for.zero( OmegaTest, out$Omega,tag= "mKrigOmega")
test.for.zero( log(det(OmegaTest)), out$lnDetOmega,
              tag="lnDetOmega")
test.for.zero( log( det( Mc)), out$lnDetCov, tag="lnDetMc" )

# check that det adjustment really works.

set.seed( 323)
x<- matrix( runif( 20), 10, 2)
temp<-  matrix( NA, 50,8)
aRangeGrid<- seq( .1,.5, ,50)
lambdaGrid<- 10**(runif( 50, -2,0))
Q<- qr.qy( qr( cbind( rep(1,10),x) ), diag( 1,10))
Q2<- Q[,4:10]
y<- rnorm(10)

testDet<- function(lambda, obj) 
{
  D2 <- obj$matrices$D[obj$matrices$D > 0]
  u2 <- obj$matrices$u[obj$matrices$D > 0]
  lD <- D2 * lambda
  N2 <- length(D2)
  sigma.MLE <- (sum((D2 * (u2)^2)/(1 + lD)))/N2
  lnDetCov <- -sum(log(D2/(1 + lD)))
#  -1 * (-N2/2 - log(2 * pi) * (N2/2) - (N2/2) * log(sigma.MLE) - 
#          (1/2) * lnDetCov)
  return( c(lnDetCov, sigma.MLE) )
}

for ( k in 1:50) {
  out<- mKrig( x,y,  aRange = aRangeGrid[k], 
                    lambda = lambdaGrid[k] 
               )
# turn off warnings for lambda search because all we want are
# matrix decompositions independent of lambda
  out2<- Krig( x,y, aRange= aRangeGrid[k],
               cov.args=list( Covariance = "Exponential"),
                            give.warnings=FALSE)
             
  Mc<-   exp( -rdist( x,x)/aRangeGrid[k] ) + lambdaGrid[k]* diag( 1,10)
  X<- out$Tmatrix
  temp[k,]<-c(
      out$lnDetCov,
      out$lnDetOmega,
      log( det( solve(t( Q2)%*%Mc%*%Q2) ) ),
      log( det(Mc) ),
      -1*log( det( t(X)%*%solve(Mc)%*%X ) ),
    testDet( lambdaGrid[k], out2 ),
      out$summary["sigma2"]
      )
}


test.for.zero( temp[,2], temp[,5], tag="testing det Omega formula")

resid<- temp[,1]  - temp[,2] + temp[,3]
test.for.zero( mean(resid), resid, relative=FALSE,
               tag="REML Det shortcut")
#### testing Krig verses mKrig
#
test.for.zero( temp[,3], -temp[,6], 
               tag="Q2 Det and Eigen Det")
###### testing sigma.MLE from mKrig and Krig

test.for.zero( (7/10)*temp[,7],  temp[,8], 
               tag="sigma.MLE Krig verses mKrig")


