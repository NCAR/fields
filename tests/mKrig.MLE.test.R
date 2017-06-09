 # fields is a package for analysis of spatial data written for
  # the R software environment .
  # Copyright (C) 2017
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


library( fields )
options( echo=FALSE)
test.for.zero.flag<- 1

#
##### generate test data
#

genCovMat = function(x, theta, lambda) {
  distanceMatrix<- rdist(x,x)
  Sigma<- exp( -distanceMatrix/theta ) + diag(x=lambda, nrow=nrow(distanceMatrix))
  return(Sigma)
}

#generate observation locations
n=500
x = matrix(runif(2*n), nrow=n)

#generate observations at the locations
trueTheta = .2
trueLambda = .1
Sigma = genCovMat(x, trueTheta, trueLambda)

U = chol(Sigma)
y = t(U)%*%as.vector(rnorm(n))

#
######set MLE computation parameters
#

testThetas = seq(from=trueTheta/2, to=2*trueTheta, length=20)
par.grid=list(theta=testThetas)
guessLambda = trueLambda

#
##### test using distance matrix
#

print("testing using distance matrix")

set.seed(1)
out1 = mKrig.MLE(x, y, lambda=guessLambda, par.grid=par.grid,
                 cov.args= list(Distance="rdist"))
lambda.MLE = out1$lambda.MLE
theta.MLE = out1$cov.args.MLE$theta

#perform mKrig at MLE parameters
out1 = mKrig(x, y, lambda=lambda.MLE, theta=theta.MLE, cov.args= list(Distance="rdist"))
print("finished default case")

set.seed(1)
out2 = mKrig.MLE(x, y, lambda=guessLambda, par.grid=par.grid)
lambda.MLE = out2$lambda.MLE
theta.MLE = out2$cov.args.MLE$theta

#perform mKrig at MLE parameters
out2 = mKrig(x, y, lambda=lambda.MLE, theta=theta.MLE)

print("finished compact distance matrix case")

#
##### test comatibility with other fields functions
#

temp1<- predict( out1)
temp2<- predict( out2)
test.for.zero( temp1, temp2, tag="predict compatibility: rdist with compact versus normal rdist")

#
##### test SE
#

temp1 = predictSE(out1)
temp2 = predictSE(out2)

test.for.zero( temp1, temp2, tag="predictSE compatibility: rdist with compact versus normal rdist")





cat("all done with mKrig.MLE tests", fill=TRUE)
options( echo=TRUE)
