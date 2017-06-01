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
#ind<- !is.na( y)
#x<- x[ind,]
#y<- y[ind]

################ test that optim results match the model evaluated 
################ at the optimized parameters.
optim.args = list(method = "BFGS", 
                  control = list(fnscale = -1, parscale = c(0.5, 0.5), 
                                 ndeps = c(0.05,0.05)))

MLEfit0 <- mKrigMLEJoint(x, y, lambda.start=.5, 
                         cov.params.start= list(theta=1.2), 
                         cov.fun="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern", smoothness=1.0),
                         na.rm=TRUE,
                       mKrig.args = list( m=1),
                         verbose=FALSE)
test.for.zero( MLEfit0$summary["lnProfileLike.FULL"], MLEfit0$optimResults$value)
 

obj0<- mKrig( x,y, cov.args = list(Covariance = "Matern",
                                smoothness = 1.0),
                                na.rm=TRUE, m=1,
                                lambda= MLEfit0$pars.MLE[1],
                                theta=MLEfit0$pars.MLE[2])
test.for.zero( MLEfit0$summary["lnProfileLike.FULL"],
                  obj0$lnProfileLike.FULL)

test.for.zero( MLEfit0$summary["rhoMLE"],obj0$rho.MLE)

par.grid<- list( theta= c(.99, 1.0, 1.01)*MLEfit0$summary["theta"] )
MLEfit1<- mKrigMLEGrid(x, y,  
                           cov.fun = "stationary.cov", 
                         cov.args  = list(Covariance = "Matern", smoothness = 1.0),
                          par.grid = par.grid, 
                            lambda = .5, 
                    lambda.profile = TRUE, 
                           mKrig.args = list( m=1),
                           na.rm=TRUE,
                           verbose = FALSE) 

hold<- (MLEfit1$summary[1,"lnProfileLike.FULL"] < MLEfit1$summary[2,"lnProfileLike.FULL"]) &
  (MLEfit1$summary[3,"lnProfileLike.FULL"] < MLEfit1$summary[2,"lnProfileLike.FULL"])

test.for.zero(as.numeric(hold), 1, relative=FALSE)



lambdaGrid<-  c(.99, 1.0, 1.01)*MLEfit0$summary["lambda"]
par.grid<- list( theta=  rep(MLEfit0$summary["theta"] ,3 ) )
MLEfit2 <- mKrigMLEGrid(x, y, 
                                 cov.fun = "stationary.cov", 
                                 cov.args  = list(Covariance = "Matern", smoothness = 1.0),
                                 mKrig.args = list( m=1),
                                 par.grid = par.grid, 
                                 lambda = lambdaGrid, 
                                 lambda.profile = FALSE, 
                                 verbose = FALSE)
hold<- (MLEfit2$summary[1,"lnProfileLike.FULL"] < MLEfit2$summary[2,"lnProfileLike.FULL"]) &
  (MLEfit2$summary[3,"lnProfileLike.FULL"] < MLEfit2$summary[2,"lnProfileLike.FULL"])
test.for.zero(as.numeric(hold), 1, relative=FALSE)

MLEfit3<- MLESpatialProcess( x,y, 
                             cov.args  = list(Covariance = "Matern", smoothness = 1.0),
                             mKrig.args = list( m=1)
                             )

test.for.zero(MLEfit0$summary[1:5], 
              (MLEfit3$MLEJoint$summary[1:5]), tol=2e-3 )

obj<- spatialProcess( x, y, mKrig.args= list(m = 1),
                         theta = MLEfit0$summary[3] )

obj1<- spatialProcess( x, y, mKrig.args= list(m = 1)
                       )

test.for.zero(MLEfit0$summary[1], 
              obj$lnProfileLike.FULL )

test.for.zero(MLEfit0$summary[1], 
              obj1$lnProfileLike.FULL)


# testing Krig function 

out1<- Krig( x,y,  cov.fun="stationary.cov",
            
            cov.args = list(Covariance = "Matern",
                            smoothness=1.0, theta=.9),
            na.rm=TRUE,
              m=2)

genCovMat = function(x, theta, lambda) {
  distanceMatrix<- rdist(x,x)
  Sigma<- Matern( distanceMatrix/theta, smoothness=1.0 ) + diag(x=lambda, nrow=nrow(distanceMatrix))
  return(Sigma)
}

#generate observation locations
set.seed( 223)
n=50
x = matrix(runif(2*n), nrow=n)
#generate observations at the locations
trueTheta = .1
trueLambda = .1
Sigma = genCovMat(x, trueTheta, trueLambda)

U = chol(Sigma)
M<- 1e4
set.seed( 332)
y = t(U)%*%matrix( rnorm(n*M), n,M)

optim.args = list(method = "BFGS", 
                  control = list(fnscale = -1, parscale = c(0.5, 0.5), 
                                 ndeps = c(0.05,0.05)))

MLEfitA <- mKrigMLEJoint(x, y, lambda.start=.5, 
                         cov.params.start= list(theta=.12), 
                         cov.fun="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern",
                                         smoothness=1.0),
                         na.rm=TRUE,
                         mKrig.args = list( m=0),
                         verbose=FALSE)
test.for.zero( MLEfitA$summary["lambda"],.1, tol=.02)
test.for.zero( MLEfitA$summary["theta"],.1, tol=.02)
test.for.zero( MLEfitA$summary["rhoMLE"], 1.0, tol=.002)

### now test REML fitting
MLEfitB <- mKrigMLEJoint(x, y, lambda.start=.5, 
                         cov.params.start= list(theta=.12), 
                         cov.fun="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern",
                                         smoothness=1.0),
                         na.rm=TRUE,
                         mKrig.args = list( m=0),
                         REML=TRUE,
                         verbose=FALSE)

test.for.zero( MLEfitB$summary["lambda"],.1, tol=.02)
test.for.zero( MLEfitB$summary["theta"],.1, tol=.02)
test.for.zero( MLEfitB$summary["rhoMLE"], 1.0, tol=.002)



MLEfitC <- mKrigMLEJoint(x, y, lambda.start=.5, 
                         cov.params.start= list(theta=.12), 
                         cov.fun="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern",
                                         smoothness=1.0),
                         na.rm=TRUE,
                         mKrig.args = list( m=2),
                         REML=FALSE,
                         verbose=FALSE
                         )
          
test.for.zero( MLEfitC$summary["lambda"],  .1, tol=.02)
test.for.zero( MLEfitC$summary[ "theta"],  .1, tol=.02)
test.for.zero( MLEfitC$summary["rhoMLE"], 1.0, tol=.002)


MLEfitA$summary
MLEfitB$summary
MLEfitC$summary


#generate observation locations
set.seed( 223)
NS<- 25
hold<- matrix(NA, nrow=NS, ncol=7 )

n=500

x = matrix(runif(2*n), nrow=n)
trueTheta = .1
trueLambda = .04
Sigma = genCovMat(x, trueTheta, trueLambda)
U = chol(Sigma)
set.seed( 332)
for( k in 1:NS){
cat(k, " ")
#generate observations at the locations
M<- 1
y = t(U)%*%matrix( rnorm(n*M), n,M)

MLEfitC <- mKrigMLEJoint(x, y, lambda.start=.5, 
                         cov.params.start= list(theta=.12), 
                         cov.fun="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern",
                                         smoothness=1.0),
                         na.rm=TRUE,
                         mKrig.args = list( m=2),
                         REML=FALSE,
                         verbose=FALSE)

hold[k,]<- MLEfitC$summary
}


cat("all done with mKrigMLEGrid tests", fill=TRUE)
options( echo=TRUE)


