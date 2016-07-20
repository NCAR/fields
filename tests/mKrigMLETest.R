# Test adapted from fields package, under GPL license

library( fields )
options( echo=FALSE)

#
##### generate test data
#

data( ozone2)
# x is a two column matrix where each row is a location in lon/lat 
# coordinates
x<- ozone2$lon.lat
# y is a vector of ozone measurements at day 16 a the locations. 
y<- ozone2$y[16,]
ind<- !is.na( y)
x<- x[ind,]
y<- y[ind]


optim.args = list(method = "BFGS", 
                  control = list(fnscale = -1, parscale = c(0.5, 0.5), 
                                 ndeps = c(0.05,0.05)))

MLEfit0 <- mKrigMLEJoint(x, y, lambda.guess=.5, 
                         cov.params.guess= list(theta=1.2), 
                         cov.fun="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern",
                                         smoothness = 1.0) ,
                         na.rm=TRUE, verbose=TRUE)
test.for.zero( MLEfit0$summary["lnProfile"], MLEfit0$optimResults$value)

par.grid<- list( theta=seq( .2, 4,,10) )
MLEfit1<- mKrigMLEGrid(x, y,  
                           cov.fun = "stationary.cov", 
                         cov.args  = list(Covariance = "Matern", smoothness = 1.0),
                          par.grid = par.grid, 
                            lambda = .5, 
                    lambda.profile = TRUE, 
                           verbose = FALSE) 
lambda.start<- 10^(seq( -2,1,,10)) 
par.grid<- list( theta=  rep(.15 , 10) )
MLEfit2 <- mKrigMLEGrid(x, y, 
                                 cov.fun = "stationary.cov", 
                                 cov.args  = list(Covariance = "Matern", smoothness = 1.0),
                                 par.grid = par.grid, 
                                 lambda = lambda.start, 
                                 lambda.profile = FALSE, 
                                 verbose = TRUE) 
MLEfit3<- MLESpatialProcess( x,y)
obj<- spatialProcess( x,y, verbose=FALSE)
test.for.zero(MLEfit0$summary[1:5], (MLEfit3$MLEInfo$MLEJoint$summary)[1:5], tol=2e-3 )




genCovMat = function(x, theta, lambda) {
  distanceMatrix<- rdist(x,x)
  Sigma<- Matern( distanceMatrix/theta, smoothness=1.0 ) + diag(x=lambda, nrow=nrow(distanceMatrix))
  return(Sigma)
}


#generate observation locations
set.seed( 223)
n=200
x = matrix(runif(2*n), nrow=n)

#generate observations at the locations
trueTheta = .1
trueLambda = .1
Sigma = genCovMat(x, trueTheta, trueLambda)

U = chol(Sigma)
M<- 1e4
set.seed( 332)
y = t(U)%*%matrix( rnorm(n*M), n,M)



MLEfit <- MLESpatialProcess(x, y, 
                            cov.function = "stationary.cov",
                            cov.args = list(Covariance = "Matern",
                                            smoothness = 1.0) , m=1)
print( MLEfit$summary)
print( MLEfit$M)

testMLE<- mKrig(x, y, 
                cov.function = "stationary.cov",
                cov.args = list(Covariance = "Matern", theta= 
                                smoothness = 1.0) , m=1 )

data( ozone2)
# x is a two column matrix where each row is a location in lon/lat 
# coordinates
x<- ozone2$lon.lat
# y is a vector of ozone measurements at day 16 a the locations. 
y<- ozone2$y[16,]
ind<- !is.na( y)
x<- x[ind,]
y<- y[ind]


optim.args = list(method = "BFGS", 
                  control = list(fnscale = -1, parscale = c(0.5, 0.5), 
                                 ndeps = c(0.05,0.05)))

MLEfit0 <- mKrigMLEJoint(x, y, lambda.guess=.5, 
                         cov.params.guess= list(theta=1.2), 
                         cov.fun="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern",
                                            smoothness = 1.0) ,
                         na.rm=TRUE, verbose=TRUE)
                           
obj<- spatialProcess( x,y, verbose=FALSE)

#
######set MLE computation parameters
#

testThetas = seq(from=trueTheta/2, to=2*trueTheta, length=20)
par.grid=list(theta=testThetas)
guessLambda = trueLambda

#
##### test using distance matrix

out1 = mKrigMLEGrid(x, y, lambda=guessLambda,
                    par.grid=list(theta=testThetas), 
                    cov.args= list(Distance="rdist"))

plot( out1$par.grid[,1], out1$summary[,2])
xg<- seq( min(testThetas), max(testThetas),,2000)
ind<- which.max( splint( out1$par.grid[,1], out1$summary[,2] , xg) )
thetaMLE<- xg[ind]
#perform mKrig at MLE parameters
out2 = mKrig(x, y, lambda=out1$lambda.MLE,
               theta=out1$cov.args.MLE$theta,
               cov.args= list(Distance="rdist"))
                   

test.for.zero( out2$lnProfileLike.FULL, max( out1$summary[,2]) )


cat("all done with mKrigMLEGrid tests", fill=TRUE)
options( echo=TRUE)


