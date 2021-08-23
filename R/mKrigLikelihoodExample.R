library( fields)
# load newer version of mKrig function
source("mKrig.R") # patch for fields 12.5
# generate an exponential spatial sample
#
n<- 500
set.seed(122)
s<- matrix( runif(n*2), n,2)
Sigma<- Exp.cov( s,s,aRange=.2)
# assume sigma2 = 1
g<- t(chol(Sigma))%*%rnorm( n)
tauTrue<- .05
y<- g + tauTrue* rnorm( n)

m1<- 30
m2<- 25
sigma2Grid<- seq( .2,1.5,length.out=m1)
aRangeGrid<- seq( .05,.4, length.out=m2)

logLike<- matrix( NA, m1, m2)
# grid search over the sigma2 and range with 
# nugget fixed (in this case cheating because
# we know the true value!)

for ( j in 1:m1){
  cat( j, " ")
  for( k in 1:m2){
    logLike[j,k]<- mKrig( s,y,cov.args= list( aRange= aRangeGrid[k]),
                     sigma2= sigma2Grid[j], tau=tauTrue)$lnLike.FULL
  }
}

llSurface<- list( x=sigma2Grid, y=aRangeGrid, z=logLike )
# handy way to find the maximum of the image object
hold<- which.max.image( llSurface)
print( hold)

# take a look at the surface
image.plot( llSurface, xlab="sigma2", ylab="range parameters")
title("log likelihood surface  exp. covariance")
contour( add=TRUE,llSurface, col="grey" )
points( 1.0, .2, col="white", pch=16, cex=1.2)
points( hold$x, hold$y, pch=16, col="grey")

# 95% confidence set  using -2 log likelihood  is approx   chi squared 
CIlevel<- hold$z - .5* qchisq(.95,2)
contour( add=TRUE,llSurface, col="grey", labels="", lwd=3,  levels= CIlevel)
###### fit the model now at these MLEs

fit<- mKrig( s, y, cov.args= list( aRange= hold$y),
             sigma2= hold$x, tau=tauTrue)
print( fit )



