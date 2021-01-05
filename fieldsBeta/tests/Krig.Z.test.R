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

#cat("A very nasty case with knots and weights",fill=TRUE)

set.seed(123)
x<- matrix( runif( 30), 15,2)
Z<- matrix( rnorm(30), 15,2)
y<- rnorm( 15)*.01 + 5*(x[,1]**3 +  (x[,2]-.5)**2) +  (Z[,1] +Z[,2])*.001
#weights<- runif(15)*10

# first without knots compare default to fixed

out.new<-  Krig( x,y,Z=Z, cov.function=Exp.cov, give.warnings=FALSE)-> out.new

out.new2<- Krig( x,y,Z=Z, cov.function=Exp.cov,
                 lambda=1)
 

##########
## compute test using linear algebra

K<- Exp.cov( x,x)
lambda<-1
M<- (lambda* diag(nrow( x)) + K)
T<- cbind( rep(1,15), x, Z)
temp.d<- c(solve( t(T) %*% solve( M)%*%T) %*% t(T)%*% solve( M) %*% y)
temp.c<- solve( M)%*% ( y - T%*% temp.d)

# test for d coefficients
test.for.zero( out.new2$d, temp.d, tag=" d coef")
# test for c coefficents
test.for.zero( out.new2$c, temp.c, tag="c coef" )


####### testing predict function 
hold2<- predict( out.new2, x=x, Z=Z, just.fixed=TRUE)
hold3<- predict( out.new2, x=x, Z=Z, drop.Z=TRUE)
hold4<- predict( out.new2, x=x, Z=Z, drop.Z=TRUE, just.fixed=TRUE)

hold<-T%*%temp.d
test.for.zero( hold, hold2, tag="predict for null fixed" )

hold<-T[,1:3]%*%temp.d[1:3] + K %*% temp.c
test.for.zero( hold, hold3, tag="predict for null spatial"  )

hold<-T[,1:3]%*%temp.d[1:3]
test.for.zero( hold, hold4, tag="predict for null drift" )

######tests where coefficients  are recomputed from object
hold2<- predict( out.new,y=y, lambda=1.0, x=x, Z=Z, just.fixed=TRUE)
hold3<- predict( out.new,y=y, lambda=1.0, x=x, Z=Z, drop.Z=TRUE)
hold4<- predict( out.new, y=y, lambda=1.0, x=x, Z=Z, 
                      drop.Z=TRUE, just.fixed=TRUE)

hold<-T%*%temp.d
test.for.zero( hold, hold2, tag="predict for null fixed" )

hold<-T[,1:3]%*%temp.d[1:3] + K %*% temp.c
test.for.zero( hold, hold3, tag="predict for null spatial" )

hold<-T[,1:3]%*%temp.d[1:3]               
test.for.zero( hold, hold4, tag="predict for null drift " )
#

####### tests using predict.se
 x<- ChicagoO3$x
 y<- ChicagoO3$y
 Zcov<-  x[,1]**3 + x[,2]**3


tps.fit<-Tps( x,y, scale.type="unscaled", Z= Zcov)

# here is lazy way to get a grid.list
   fields.x.to.grid( x, nx=20,ny=20)-> gridlist

   xg<- make.surface.grid(gridlist)
   Zcov.grid<- xg[,1]**3 + xg[,2]**3

########### tests on just predict have been commented out to 
########### indicate that they are redundant given 
########### previous tests however, they could be useful for 
########### future debugging ...

# full surface with covariate
#   curv.mean1 <- predictSurface(tps.fit, grid.list = gridlist, extrap = TRUE,
##        Z =Zcov.grid, drop.Z = FALSE)$z

# just the spline surface
#   curv.mean2 <- predictSurface(tps.fit, grid.list = gridlist,
#                   extrap = TRUE,drop.Z=TRUE)$z

# explicitly here is the difference surface of curv.mean1 and curv.mean2
#   curv.mean0<- as.surface( gridlist, Zcov.grid* tps.fit$d[4])$z
#   test.for.zero( curv.mean1- curv.mean2, curv.mean0)

## new tests

   predictSurfaceSE( tps.fit, grid.list=gridlist, extrap=TRUE,
               drop.Z=TRUE)$z-> curv.var1

   predictSE( tps.fit, xg, drop.Z=TRUE)-> curv.var2
     test.for.zero( curv.var1, curv.var2)

# SE with covariates included
   predictSE( tps.fit, xg, Z=Zcov.grid, drop.Z=FALSE)**2-> curv.var1
#   as.surface( gridlist, curv.var1)$z-> curv.var1

# SE for just the spline part
   predictSE( tps.fit, xg, drop.Z=TRUE)**2-> curv.var2
#   as.surface( gridlist, curv.var2)$z-> curv.var2

# SE for just the fixed part
 predictSE( tps.fit, xg,Z=Zcov.grid, drop.Z=FALSE,
                            just.fixed=TRUE )**2-> curv.var3
# as.surface( gridlist, curv.var3)$z-> curv.var3


# calculating from more basic functions
## these tests assume that Krig.Amatrix is working correctly!

out<- tps.fit

A<- Krig.Amatrix( tps.fit,x= xg, drop.Z=TRUE)
Sigma<- out$sigmahat*Rad.cov( out$x, out$x, p=2)
S0<- out$sigmahat*Rad.cov(xg, xg, p=2)
S1<- out$sigmahat*Rad.cov(out$x, xg, p=2)

#yhat= Ay
#var( f0 - yhat)=    var( f0) - 2 cov( f0,yhat)+  cov( yhat)

   look<- S0 - t(S1)%*% t(A) - A%*%S1 +
       A%*% ( Sigma + diag(out$tauHat.MLE**2/out$weightsM))%*% t(A)
   look<- diag( look)
   test.for.zero(curv.var2 ,look,tag="SE w/o covariate")


A<- Krig.Amatrix( tps.fit,x= xg, drop.Z=FALSE,Z=Zcov.grid)
# see tps.fit$args for value of p
Sigma<- out$sigmahat*Rad.cov( out$x, out$x, p=2)
S0<- out$sigmahat*Rad.cov(xg, xg, p=2)
S1<- out$sigmahat*Rad.cov(out$x, xg, p=2)

#yhat= Ay
#var( f0 - yhat)=    var( f0) - 2 cov( f0,yhat)+  cov( yhat)

   look<- S0 - t(S1)%*% t(A) - A%*%S1 +
       A%*% ( Sigma + diag(out$tauHat.MLE**2/out$weightsM))%*% t(A)
   look<- diag( look)
   test.for.zero(curv.var1 ,look,tag="SE with covariate")


A<- Krig.Amatrix( tps.fit,x= xg, drop.Z=FALSE,Z=Zcov.grid, just.fixed=TRUE)
# see tps.fit$args for value of p
Sigma<- out$sigmahat*Rad.cov( out$x, out$x, p=2)
S0<- out$sigmahat*Rad.cov(xg, xg, p=2)
S1<- out$sigmahat*Rad.cov(out$x, xg, p=2)

#yhat= Ay
#var( f0 - yhat)=    var( f0) - 2 cov( f0,yhat)+  cov( yhat)

   look<- S0 - t(S1)%*% t(A) - A%*%S1 +
       A%*% ( Sigma + diag(out$tauHat.MLE**2/out$weightsM))%*% t(A)
   look<- diag( look)
   test.for.zero(curv.var3 ,look, tag="SE for fixed part")

cat("All done with Z tests and Krig/Tps including predict and predictSE !", 
          fill=TRUE)
options( echo=TRUE)
