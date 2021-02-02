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
#
#
#  test of fixed lambda case
#  Check against linear algebra
#

options( echo=FALSE)
test.for.zero.flag<-1

fit<- Krig( ChicagoO3$x, ChicagoO3$y, aRange=50)

x<- ChicagoO3$x
K<- Exp.cov(x, x,aRange=50)
T<- fields.mkpoly(x, 2)
W<- diag( 20)
 lambda<- fit$lambda
M<- (lambda* diag(20) + K) 
###########################
test.d<- c(solve( t(T) %*% solve( M)%*%T) %*% t(T)%*% solve( M) %*% fit$yM)
test.c<- solve( M)%*% ( fit$yM - T%*% test.d)

#compare to  fit$d
test.for.zero( test.d, fit$d, tag="Compare d coef" )
#compare to  fit$d
test.for.zero( test.c, fit$c,tag="Compare c coef" )

fit2<- Krig( ChicagoO3$x, ChicagoO3$y, aRange=50, lambda= fit$lambda)
#compare to  fit$d
test.for.zero( test.d, fit2$d, tag="Compare d coef fixed lambda" )
#compare to  fit$d
test.for.zero( test.c, fit2$c,tag="Compare c coef fixed lambda" )

# test of Krig.coef

Krig.coef( fit)->test
test.for.zero( test.d, test$d, tag="d coef Krig.coef" )
test.for.zero( test.c, test$c, tag= "c coef Krig.coef" )

Krig.coef( fit2)->test
test.for.zero( test.d, test$d,tag="d coef Krig.coef fixed" )
test.for.zero( test.c, test$c, tag="c coef Krig.coef fixed" )
# checking A matrix in the case of noreps

set.seed( 222)
weights<-  10+ runif( length(ChicagoO3$y))
#weights<- rep( 1, 20)
test2<- Krig( ChicagoO3$x, ChicagoO3$y, aRange=50, weights= weights)
Atest<- Krig.Amatrix( test2)
K<-Exp.cov(ChicagoO3$x, ChicagoO3$x,aRange=50)
H<- matrix(0, 23,23)
H[(1:20)+3 , (1:20)+3]<- K
X<- cbind( fields.mkpoly( ChicagoO3$x, 2), K)
lambda<- test2$lambda
 Alam <-  X%*%solve(
                 t(X)%*%diag(weights)%*%X + lambda*H
                 )%*% t(X)%*%diag(weights) 
 test.for.zero( Alam, Atest, tag="Amatrix no reps", tol=5e-8)

# test for new y fixed case
set.seed( 123)
ynew<- rnorm( fit2$N)

test.d<- c(solve( t(T) %*% solve( M)%*%T) %*% t(T)%*% solve( M) %*% ynew)
test.c<- solve( M)%*% ( ynew - T%*% test.d)

test<- Krig.coef( fit, y= ynew)
test.for.zero( test.d, test$d, tag= "d coef new y" )
test.for.zero( test.c, test$c, tag="c coef new y" )


Krig.coef( fit2, y= ynew)->test
test.for.zero( test.d, test$d, tag= "d coef new y fixed" )
test.for.zero( test.c, test$c, tag=" c coef new y fixed"  )

# test for multiple new y's
Krig.coef( fit2, y= cbind( ynew+ rnorm(fit2$N), ynew))->test2
test.for.zero( test.d, test2$d[,2], tag= "d coef several new y fixed" )
test.for.zero( test.c, test2$c[,2], tag=" c coef several new y fixed"  )


#cat("done with simple Krig data", fill=TRUE)


# These tests are about whether decompositions 
# handle just a fixed lambda or are more general 

# checking passing lambda or df to Krig

out<- Tps( ChicagoO3$x, ChicagoO3$y,lambda=.001 )
out2<- predict( out, lambda=.001)
test.for.zero( out2, predict( out), tag="Tps with fixed lam")

out<- Tps( ChicagoO3$x, ChicagoO3$y, df=5)
out2<- predict( out, df=5)
test.for.zero( out2, predict( out), tag="Tps with fixed df")

# same for Krig

out0<- Krig( ChicagoO3$x, ChicagoO3$y, aRange=50,lambda=.5)
out<- Krig( ChicagoO3$x, ChicagoO3$y, aRange=50,lambda=.5,GCV=TRUE)
test.for.zero( 
      predict(out0), predict( out), tag="Krig with fixed lam argument")

#A very nasty case with knots and weights

set.seed(123)
x<- matrix( runif( 30), 15,2)
y<- rnorm( 15)*.01 + x[,1]**2 +  x[,2]**2

weights<- runif(15)*10

# compare to 
Krig( x,y, cov.function=Exp.cov, weights=weights)-> out.new
Krig( x,y, cov.function=Exp.cov, weights=weights, 
          lambda=1)-> out.new2

# compute test using linear algebra
K<- Exp.cov( x, x)
H<- matrix(0, 18,18)
H[4:18, 4:18]<- K
X<- cbind( fields.mkpoly( x, 2), Exp.cov( x, x))
lambda<-1


c(   solve(t(X)%*%(weights*X) + lambda*H)%*% t(X)%*% (weights*y) )-> temp
temp.c<- temp[4:18]
temp.d<- temp[1:3]


# test for d coefficients
test.for.zero( out.new2$d, temp.d, tag=" d coef")
# test for c coefficents
test.for.zero( out.new2$c, temp.c, tag="c coef" )

# and 
test<- Krig.coef( out.new2, lambda=1)

# test for d coefficients
test.for.zero( temp.d, test$d, tag= "d fixed case")
# test for c coefficents 
test.for.zero( temp.c, test$c, tag=" c fixed case" )


ynew<- 1:15

#compare 
test<- Krig.coef( out.new, lambda=.5, y=ynew) 
test2<- Krig( x,ynew, cov.function=Exp.cov,
              lambda= .5, weights=weights)
# test for d coefficients
test.for.zero(  test2$d,test$d, tag=" d new y")
# test for c coefficents 
test.for.zero( test2$c, test$c,tag= "c new y" )




#cat("test with reps" , fill=TRUE)
#

##################################
#cat( "test  A matrix",fill=TRUE)
##################################
 
set.seed(133)
x<- matrix( runif( 30), 15,2)*2  
y<- rnorm( nrow( x))*.5 + + x[,1]**2 +  x[,2]**2
# perturb so that this example does not generate (harmless) warnings in gcv search
weights<- runif( nrow( x))*10

out.new<- Krig( x,y, weights= weights)


testY<- predict( out.new)
testY2<-  Krig.Amatrix(out.new)%*% y
test.for.zero( testY, testY2, tag="testing A matrix")

set.seed(333)
yNew<- rnorm( 15)
testY<- predict( out.new, y=yNew)
testY2<-  Krig.Amatrix(out.new)%*% yNew
test.for.zero( testY, testY2, tag="testing A matrix new y")

Alam<- Krig.Amatrix(out.new)
trA<- sum( diag( Alam))
test.for.zero( trA, out.new$eff.df,  tag="checking trace")

####### checking GCV
#######

MSE<- mean( out.new$residuals^2 * weights)
n<- length( y)
GCV<- MSE/ (1- trA/n )^2
 
test.for.zero(GCV, out.new$lambda.est[6,3], tag="GCV" )
 




