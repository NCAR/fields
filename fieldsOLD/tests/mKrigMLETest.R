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

suppressMessages(library( fields ))
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

x<-x[1:31,]
y<-y[1:31]
y[31]<-NA


################ test that optim results match the model evaluated 
################ at the optimized parameters.
optim.args = list(method = "BFGS", 
                  control = list(fnscale = -1, parscale = c(0.5, 0.5), 
                                 ndeps = c(0.05,0.05)))

MLEfit0 <- mKrigMLEJoint(x, y,  
                         cov.params.start= list(lambda=.5, theta=1.2), 
                         cov.fun="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern", smoothness=1.0),
                         na.rm=TRUE,
                       mKrig.args = list( m=1),
                         verbose=FALSE)
# check agreement with a fast return  note cov.params.start and lambda.fixed 
# are the switches to indicate this case

MLEfit0C <- mKrigMLEJoint(x, y,
                          cov.params.start = NULL,  
                                   cov.fun = "stationary.cov",
                                  cov.args = list(Covariance = "Matern",
                                                  lambda = MLEfit0$pars.MLE["lambda"],
                                                   smoothness = 1.0,
                                                   theta = MLEfit0$pars.MLE["theta"]),
                                     na.rm = TRUE,
                                mKrig.args = list( m=1),
                                   verbose = FALSE
                          )



test.for.zero( MLEfit0$summary["lnProfileLike.FULL"], MLEfit0C$summary["lnProfileLike.FULL"],
               tag="Likelihood Values optim and the fast return")
 
test.for.zero( MLEfit0$summary["lnProfileLike.FULL"], MLEfit0$optimResults$value,
                tag="Likelihood Values summary and optim")
 

obj0<- mKrig( x,y, cov.args = list(Covariance = "Matern",
                                smoothness = 1.0),
                                na.rm=TRUE, m=1,
                                lambda= MLEfit0$pars.MLE["lambda"],
                                theta=MLEfit0$pars.MLE["theta"])
test.for.zero( MLEfit0$summary["lnProfileLike.FULL"],
                  obj0$lnProfileLike.FULL,
               tag="Likelihood Values summary and direct mKrig call")

test.for.zero( MLEfit0$summary["rhoMLE"],obj0$rho.MLE, 
               tag="... and rho.MLE")

# test that grid seraching is correct
theta.MLE<- MLEfit0$summary["theta"]
par.grid<- list( theta= c(.5, 1.0, 1.5)*theta.MLE )
MLEfit1<- mKrigMLEGrid(x, y,  
                           cov.fun = "stationary.cov", 
                         cov.args  = list(Covariance = "Matern",
                                          smoothness = 1.0
                                          ),
                          par.grid = par.grid, 
                        mKrig.args = list( m=1),
                             na.rm = TRUE,
                           verbose = FALSE,
                  cov.params.start = list( lambda = .2)
                  ) 


hold<- (MLEfit1$summary[1,"lnProfileLike.FULL"] < MLEfit1$summary[2,"lnProfileLike.FULL"]) &
  (MLEfit1$summary[3,"lnProfileLike.FULL"] < MLEfit1$summary[2,"lnProfileLike.FULL"])

test.for.zero(as.numeric(hold), 1, relative=FALSE,
              tag="consistency of Likelihood values")

##########################
### now evaluate on the "grid" of lambdas found by profiling
lambda.MLEs<- MLEfit1$summary[,"lambda"]
par.grid<- list( lambda = lambda.MLEs, 
                  theta = c(.5, 1.0, 1.5)*theta.MLE )
MLEfit1B<- mKrigMLEGrid(x, y,  
                         cov.fun = "stationary.cov", 
                       cov.args  = list(Covariance = "Matern", smoothness = 1.0),
                        par.grid = par.grid, 
                      mKrig.args = list( m=1),
                           na.rm = TRUE,
                         verbose = FALSE) 
tempCol<- c( "lnProfileLike.FULL",
            "lambda", "sigmaMLE","rhoMLE")
test.for.zero( as.matrix(MLEfit1$summary[,tempCol]),
               as.matrix(MLEfit1B$summary[,tempCol]),
               tag="grid search with and w/o profile")

par.grid<- list( lambda = c(.999, 1.0, 1.001)*MLEfit0$summary["lambda"],
                 theta  =  rep(MLEfit0$summary["theta"] ,3 ) )
MLEfit2 <- mKrigMLEGrid(x, y, 
                                 cov.fun = "stationary.cov", 
                                 cov.args  = list(Covariance = "Matern",
                                                  smoothness = 1.0),
                                 mKrig.args = list( m=1),
                                 par.grid = par.grid, 
                                 verbose = FALSE)
hold<- (MLEfit2$summary[1,"lnProfileLike.FULL"] < MLEfit2$summary[2,"lnProfileLike.FULL"]) &
  (MLEfit2$summary[3,"lnProfileLike.FULL"] < MLEfit2$summary[2,"lnProfileLike.FULL"])
test.for.zero(as.numeric(hold), 1, relative=FALSE, tag="crude test of maxmimum")

MLEfit3<- MLESpatialProcess( x,y, 
                             cov.args  = list(Covariance = "Matern",
                                              smoothness = 1.0),
                             mKrig.args = list( m=1),
                             cov.params.start = list( lambda =.2, theta = NA)
                             )

test.for.zero(MLEfit0$summary[1:5]/ 
              (MLEfit3$MLEJoint$summary[1:5]), 1, tol=4e-3,
               tag="Testing MLESpatialProcess ")

######### making sure spatialProcess uses parameter information correctly

obj<- spatialProcess( x, y, mKrig.args= list(m = 1),
                          theta = MLEfit0$summary["theta"] 
                      )

obj1<- spatialProcess( x, y, mKrig.args= list(m = 1), 
                       )

test.for.zero(MLEfit0$summary[1], 
              obj$lnProfileLike.FULL,
              tag="spatialProcess finding MLE " )

test.for.zero(MLEfit0$summary[1], 
              obj1$lnProfileLike.FULL, tol=5e-8,
              tag="spatialProcess given MLE " 
              )


# testing Krig function 
#out1<- Krig( x,y,  cov.fun="stationary.cov",
#            cov.args = list(Covariance = "Matern",
#                            smoothness=1.0, theta=.9),
#           na.rm=TRUE,
#              m=2)


genCovMat = function(x, theta, lambda) {
  distanceMatrix<- rdist(x,x)
  Sigma<- Matern( distanceMatrix/theta, smoothness=1.0 ) + diag(x=lambda, nrow=nrow(distanceMatrix))
  return(Sigma)
}

#generate observation locations
set.seed( 223)
n=25
x = matrix(runif(2*n), nrow=n)
#generate observations at the locations
trueTheta = .1
trueLambda = .1
Sigma = genCovMat(x, trueTheta, trueLambda)

U = chol(Sigma)
M<- 1e4 # lots of replicated fields.
set.seed( 332)
y = t(U)%*%matrix( rnorm(n*M), n,M)

out<- mKrig( x,y, lambda=trueLambda, theta=trueTheta*.1, 
       cov.function ="stationary.cov",cov.args = list(Covariance = "Matern",
                                                      smoothness=1.0)
)
       
optim.args = list(method = "BFGS", 
                  control = list(fnscale = -1, parscale = c(0.5, 0.5), 
                                 ndeps = c(0.05,0.05)))

MLEfitA <- mKrigMLEJoint(x, y, 
                         cov.params.start= list(theta=.12, lambda=.5), 
                         cov.fun="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern",
                                         smoothness=1.0),
                         na.rm=TRUE,
                         mKrig.args = list( m=0),
                         verbose=FALSE)

cat("Testing mKrigMLEJoint against true values",
    fill=TRUE)
test.for.zero( MLEfitA$summary["lambda"],.1, tol=.02)
test.for.zero( MLEfitA$summary["theta"],.1, tol=.02)
test.for.zero( MLEfitA$summary["rhoMLE"], 1.0, tol=.002)

### now test REML fitting
MLEfitB <- mKrigMLEJoint(x, y, 
                         cov.params.start= list(theta=.12, lambda=.5), 
                         cov.fun="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern",
                                         smoothness=1.0),
                         na.rm=TRUE,
                         mKrig.args = list( m=0),
                         REML=TRUE,
                         verbose=FALSE)

cat("Testing mKrigMLEJoint  with REML against true values",
    fill=TRUE)
test.for.zero( MLEfitB$summary["lambda"],.1, tol=.02)
test.for.zero( MLEfitB$summary["theta"],.1, tol=.02)
test.for.zero( MLEfitB$summary["rhoMLE"], 1.0, tol=.002)


cat("Testing mKrigMLEJoint  with REML FALSE  against true values",
    fill=TRUE)

MLEfitC <- mKrigMLEJoint(x, y,  
                         cov.params.start= list(theta=.12, lambda=.5), 
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


signif( MLEfitA$summary,5)
signif( MLEfitB$summary,5)
signif( MLEfitC$summary,5)


# simple Monte Carlo test
 
NS<- 2
n<-75  
M<- 1000


set.seed(123)
x = matrix(runif(2*n), nrow=n)
trueTheta = .1
trueLambda = .04
Sigma = genCovMat(x, trueTheta, trueLambda)
U = chol(Sigma)
set.seed( 332)
hold<- NULL
for( k in 1:NS){
cat(k, " ")
#generate observations at the locations

y = t(U)%*%matrix( rnorm(n*M), n,M)

MLEfitC <- mKrigMLEJoint(x, y,
                         cov.params.start= list(theta=.12, lambda=.5), 
                         cov.fun="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern",
                                         smoothness=1.0),
                         na.rm=TRUE,
                         mKrig.args = list( m=2),
                         REML=FALSE,
                         verbose=FALSE)

hold<- rbind( hold, c(MLEfitC$summary) )
}
cat(" ", fill=TRUE)

test.for.zero( trueTheta, mean(hold[,"theta"]), tol=6e-3,
               tag="Monte Carlo theta")
test.for.zero( trueLambda, mean(hold[,"lambda"]), tol=5e-2,
               tag="Monte Carlo lambda")


cat("all done with mKrigMLEGrid tests", fill=TRUE)
options( echo=TRUE)


