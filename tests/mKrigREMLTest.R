
library( fields)
set.seed( 123)
x<- matrix( runif( 20),10, 2)
y<- rnorm(10)


lambda<- .1
theta<- .2
out<- mKrig( x,y, theta= theta, lambda=lambda)

test.for.zero( out$lnDetOmega,
2*log( prod(diag(chol(out$Omega))))
)
    
Mc<-   exp( -rdist( x,x)/theta) + lambda* diag( 1,10)
OmegaTest<- solve(t(out$Tmatrix)%*%solve( Mc)%*% out$Tmatrix)

test.for.zero( OmegaTest, out$Omega,tag= "mKrigOmega")
test.for.zero( log(det(OmegaTest)), out$lnDetOmega,
              tag="lnDetOmega")
test.for.zero( log( det( Mc)), out$lnDetCov, tag="lnDetMc" )

# check that det adjustment really works.

set.seed( 323)
x<- matrix( runif( 20), 10, 2)
temp<-  matrix( NA, 50,8)
thetaGrid<- seq( .1,.5, ,50)
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
  rho.MLE <- (sum((D2 * (u2)^2)/(1 + lD)))/N2
  lnDetCov <- -sum(log(D2/(1 + lD)))
#  -1 * (-N2/2 - log(2 * pi) * (N2/2) - (N2/2) * log(rho.MLE) - 
#          (1/2) * lnDetCov)
  return( c(lnDetCov, rho.MLE) )
}

for ( k in 1:50) {
  out<- mKrig( x,y,  theta = thetaGrid[k], 
                    lambda = lambdaGrid[k] 
               )
  out2<- Krig( x,y, theta= thetaGrid[k],
               cov.args=list( Covariance = "Exponential" )
             )
  Mc<-   exp( -rdist( x,x)/thetaGrid[k] ) + lambdaGrid[k]* diag( 1,10)
  X<- out$Tmatrix
  temp[k,]<-c(
      out$lnDetCov,
      out$lnDetOmega,
      log( det( solve(t( Q2)%*%Mc%*%Q2) ) ),
      log( det(Mc) ),
      -1*log( det( t(X)%*%solve(Mc)%*%X ) ),
      testDet( lambdaGrid[k], out2 ),
      out$rho.MLE
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
###### testing rho.MLE from mKrig and Krig

test.for.zero( (7/10)*temp[,7],  temp[,8], 
               tag="rho.MLE Krig verses mKrig")


#lm.out<-lm( temp[,1]~ temp[,c(2:3)])
#summary( lm.out)

