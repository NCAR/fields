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
                         cov.params.start= list(lambda=.5, aRange=1.2), 
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
                              cov.function = "stationary.cov",
                                  cov.args = list(Covariance = "Matern",
                                                  lambda = MLEfit0$pars.MLE["lambda"],
                                                   smoothness = 1.0,
                                                   aRange = MLEfit0$pars.MLE["aRange"]),
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
                                aRange=MLEfit0$pars.MLE["aRange"])
test.for.zero( MLEfit0$summary["lnProfileLike.FULL"],
                  obj0$summary["lnProfileLike.FULL"],
               tag="Likelihood Values summary and direct mKrig call")

test.for.zero( MLEfit0$summary["sigma2"],obj0$summary["sigma2"], 
               tag="... and sigma^2.MLE")

# test that grid seraching is correct
aRange.MLE<- MLEfit0$summary["aRange"]
par.grid<- list( aRange= c(.5, 1.0, 1.5)*aRange.MLE )
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
                  aRange = c(.5, 1.0, 1.5)*aRange.MLE )
MLEfit1B<- mKrigMLEGrid(x, y,  
                    cov.function = "stationary.cov", 
                       cov.args  = list(Covariance = "Matern", smoothness = 1.0),
                        par.grid = par.grid, 
                      mKrig.args = list( m=1),
                           na.rm = TRUE,
                         verbose = FALSE) 
tempCol<- c( "lnProfileLike.FULL",
            "lambda", "tau","sigma2")
test.for.zero( as.matrix(MLEfit1$summary[,tempCol]),
               as.matrix(MLEfit1B$summary[,tempCol]),
               tag="grid search with and w/o profile")

par.grid<- list( lambda = c(.999, 1.0, 1.001)*MLEfit0$summary["lambda"],
                 aRange  =  rep(MLEfit0$summary["aRange"] ,3 ) )
MLEfit2 <- mKrigMLEGrid(x, y, 
                                 cov.function  = "stationary.cov", 
                                 cov.args  = list(Covariance = "Matern",
                                                  smoothness = 1.0),
                                 mKrig.args = list( m=1),
                                 par.grid = par.grid, 
                                 verbose = FALSE)
hold<- (MLEfit2$summary[1,"lnProfileLike.FULL"] < MLEfit2$summary[2,"lnProfileLike.FULL"]) &
  (MLEfit2$summary[3,"lnProfileLike.FULL"] < MLEfit2$summary[2,"lnProfileLike.FULL"])
test.for.zero(as.numeric(hold), 1, relative=FALSE, tag="crude test of maxmimum")

#MLEfit3<- MLESpatialProcess( x,y, 
#                             cov.args  = list(Covariance = "Matern",
 #                                             smoothness = 1.0),
 #                            mKrig.args = list( m=1),
 #                            cov.params.start = list( lambda =.2, aRange = NA)
 #                            )

MLEfit3<- spatialProcess( x,y, 
                          cov.args  = list(Covariance = "Matern",
                                               smoothness = 1.0),
                         mKrig.args = list( m=1),
                   cov.params.start = list( lambda =.2)
                                                      )

test.for.zero(MLEfit0$summary[1:5]/ 
              (MLEfit3$summary[1:5]), 1, tol=2e-3,
               tag="Testing MLE from spatialProcess ")

######### making sure spatialProcess uses parameter information correctly

obj<- spatialProcess( x, y, mKrig.args= list(m = 1),
                          lambda= MLEfit0$summary["lambda"], 
                          aRange = MLEfit0$summary["aRange"] 
                      )

obj1<- spatialProcess( x, y, mKrig.args= list(m = 1), 
                       )

test.for.zero(MLEfit0$summary[1], 
              obj$summary["lnProfileLike.FULL"],
              tag="spatialProcess finding MLE " )

test.for.zero(MLEfit0$summary[1], 
              obj1$summary["lnProfileLike.FULL"], tol=5e-8,
              tag="spatialProcess given MLE " 
              )


# testing Krig function 
#out1<- Krig( x,y,  cov.fun="stationary.cov",
#            cov.args = list(Covariance = "Matern",
#                            smoothness=1.0, aRange=.9),
#           na.rm=TRUE,
#              m=2)


genCovMat = function(x, aRange, lambda) {
  distanceMatrix<- rdist(x,x)
  Sigma<- Matern( distanceMatrix/aRange, smoothness=1.0 )
  + diag(x=lambda, nrow=nrow(distanceMatrix))
  return(Sigma)
}

#generate observation locations
set.seed( 22)
n=100
x = matrix(runif(2*n), nrow=n)
#generate observations at the locations
trueARange = .1
trueLambda = .1

distanceMatrix<- rdist(x,x)
Sigma<- Matern( distanceMatrix/trueARange, smoothness=1.0 )
U = chol(Sigma)
M<- 1e5 # lots of replicated fields.
set.seed( 332)
y = t(U)%*%matrix( rnorm(n*M), n,M) + 
             sqrt(trueLambda)*matrix( rnorm(n*M), n,M)

out<- mKrig( x,y, lambda=trueLambda, aRange=trueARange, 
       cov.function ="stationary.cov",cov.args = list(Covariance = "Matern",
                                                      smoothness=1.0)
)
       
optim.args = list(method = "BFGS", 
                  control = list(fnscale = -1, parscale = c(0.5, 0.5), 
                                 ndeps = c(0.05,0.05)))

MLEfitA <- mKrigMLEJoint(x, y, 
                         cov.params.start= list(aRange=.2, lambda=.01), 
                         cov.function="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern",
                                         smoothness=1.0),
                         na.rm=TRUE,
                         abstol = 1e-7,
                         mKrig.args = list( m=0),
                         verbose=FALSE)

cat("Testing mKrigMLEJoint against true values",
    fill=TRUE)
test.for.zero( MLEfitA$summary["lambda"],.1, tol=.002)
test.for.zero( MLEfitA$summary["aRange"],.1, tol=.0005)
test.for.zero( MLEfitA$summary["sigma2"], 1.0, tol=.0005)

### now test REML fitting
MLEfitB <- mKrigMLEJoint(x, y, 
                         cov.params.start= list(aRange=.12, lambda=.5), 
                         cov.function="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern",
                                         smoothness=1.0),
                         na.rm=TRUE,
                         mKrig.args = list( m=0),
                         REML=TRUE,
                         verbose=FALSE)

cat("Testing mKrigMLEJoint  with REML against true values",
    fill=TRUE)
test.for.zero( MLEfitB$summary["lambda"],.1, tol=.002)
test.for.zero( MLEfitB$summary["aRange"],.1, tol=.002)
test.for.zero( MLEfitB$summary["sigma2"], 1.0, tol=.001)


cat("Testing mKrigMLEJoint  with REML FALSE  against true values",
    fill=TRUE)

MLEfitC <- mKrigMLEJoint(x, y,  
                         cov.params.start= list(aRange=.12, lambda=.5), 
                         cov.function ="stationary.cov",
                         optim.args=optim.args,
                         cov.args = list(Covariance = "Matern",
                                         smoothness=1.0),
                         na.rm=TRUE,
                         mKrig.args = list( m=2),
                         REML=FALSE,
                         verbose=FALSE
                         )
          
test.for.zero( MLEfitC$summary["lambda"],  .1, tol=.02)
test.for.zero( MLEfitC$summary[ "aRange"],  .1, tol=.02)
test.for.zero( MLEfitC$summary["sigma2"], 1.0, tol=.002)


cat("all done with mKrigMLEGrid tests", fill=TRUE)
options( echo=TRUE)


