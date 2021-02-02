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



# this is a test script to verify the likelihood computations are 
# correct with the eigen decomposition format used in Krig
# see Krig.flplike for the concise computation.
#  

suppressMessages(library(fields))

options( echo=FALSE)
test.for.zero.flag<- 1

# utility function foor testing 
REML.test <- function(x, y, rho, sigma2, theta, nu = 1.5) {
  Tmatrix <- fields.mkpoly(x, 2)
  qr.T <- qr(Tmatrix)
  N <- length(y)
  Q2 <- qr.yq2(qr.T, diag(1, N))
  ys <- t(Q2) %*% y
  N2 <- length(ys)
  A <- (rho * Matern(rdist(x, x), range = theta, smoothness = nu) + 
          sigma2 * diag(1, N))
  A <- t(Q2) %*% A %*% Q2
  Ac <- chol(A)
  w <- backsolve(Ac, ys, transpose = TRUE)
  REML.like <- (N2/2) * log(2 * pi) + (1/2) * 2 * sum(log(diag(Ac))) + 
    (1/2) * t(w) %*% w
  REML.like <- -1 * REML.like
  ccoef <- rho * Q2 %*% solve(A) %*% ys
  return(list(REML.like = REML.like, A = A, ccoef = ccoef, 
              quad.form = t(w) %*% w, rhohat = (t(w) %*% w/N2) * rho, 
              det = 2 * sum(log(diag(Ac))), N2 = N2))
}

data( ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]
is.good <- !is.na( y)
x<- x[is.good,]
y<- y[is.good]

aRange<- 2.0

# check log likelihood calculation
  nu<- 1.5
  lambda<- .2
  out<- mKrig( x,y, aRange=aRange,Covariance="Matern", smoothness=nu, lambda=lambda) 

# peg sigma and tau as MLEs from mKrig
  sigma <- out$summary["sigma2"]
  tau2<- sigma*lambda
  N<- length( y)
  dd<- rdist( x,x)
  M<-  sigma* Matern( dd, range= aRange, smoothness=nu) + tau2* diag( 1, N)
  X<- fields.mkpoly( x, 2)
  Mi<- solve( M)
  betahat<-  solve(t(X)%*%Mi%*%X)%*% t(X)%*% Mi%*% y
  res<- y - X%*%betahat
  ccoef<- ( Mi%*% ( res))*sigma

# sanity check that estimates are the same
  test.for.zero( ccoef, out$c.coef, tag="check ccoef")

# find full log likelihood
  chol(M)-> cM
  lLike<-  -(N/2)*log(2*pi) - (1/2)* (2*sum( log( diag(cM)))) - (1/2)* t(res)%*% Mi %*% res
  test.for.zero( lLike, out$summary["lnProfileLike.FULL"], tag="llike profile from mKrig")
  
# formula for full likelihood using peices from mKrig
# lLike.test<- -(N/2)*log(2*pi) - (1/2)* out$lnDetCov - (1/2)*(N)*log( sigma) - (1/2)*out$quad.form/sigma
#  test.for.zero( lLike, lLike.test, tag="llike full verses sigmahat")
  
# REML check
  nu<- 1.5
  aRange<- .6
  obj<- Krig( x,y, aRange=aRange,Covariance="Matern", smoothness=nu )

# sanity check that c coefficients agree with Krig
  sigma<- 500
  lambda<- .2
  tau2<- lambda*sigma

  hold<- REML.test( x,y,sigma, tau2, aRange, nu=1.5)
  ccoef2<- Krig.coef( obj, lambda)$c
  test.for.zero( hold$ccoef, ccoef2, tag="ccoefs")

# check RSS with Krig decomposition.
  RSS1<- sum( (lambda*ccoef2)**2)
  lD <- obj$matrices$D * lambda
  RSS2 <- sum(((obj$matrices$u * lD)/(1 + lD))^2) 
  test.for.zero( RSS2, RSS1, tag=" RSS using matrices")

# check quadratic form with Krig
  D.temp<- obj$matrices$D[  obj$matrices$D>0]
  A3test<- (1/lambda)* obj$matrices$V %*% diag((D.temp*lambda)/ (1 +D.temp*lambda) )%*% t( obj$matrices$V)
  test.for.zero(solve(A3test), hold$A/sigma, tol=5e-8)
  Quad3<-   sum( D.temp*(obj$matrices$u[obj$matrices$D>0])^2/(1+lambda*D.temp))

  test.for.zero( hold$quad.form, Quad3/sigma, tag="quad form")

# test determinants
  N2<- length( D.temp)
  det4<- -sum( log(D.temp/(1 + D.temp*lambda)) )
  det1<- sum( log(eigen( hold$A/sigma)$values)) 
  test.for.zero( det1, det4, tag="det" )

# test REML Likelihood
  lLikeREML.test<--1*( (N2/2)*log(2*pi) - (1/2)*(sum( log(D.temp/(1 + D.temp*lambda)) ) - N2*log(sigma)) +
                                      (1/2)*sum( lD*(obj$matrices$u)^2/(1+lD)) /(lambda*sigma) )

test.for.zero( hold$REML.like, lLikeREML.test, tag="REML using matrices")




cat("all done with likelihood  tests", fill=TRUE)
options( echo=TRUE)
      
