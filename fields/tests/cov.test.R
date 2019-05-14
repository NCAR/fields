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
options( echo=FALSE)
test.for.zero.flag<- 1
data(ozone2)
y<- ozone2$y[16,]
x<- ozone2$lon.lat
#
# Omit the NAs
good<- !is.na( y)
x<- x[good,]
y<- y[good]
x1<- x[1:5,]
x2<- x[6:11,]

look<- exp(-1*rdist(x1,x2)/4)
look2<- stationary.cov( x1,x2, theta=4)
look3<- Exp.cov( x1, x2, theta=4.0)
test.for.zero( look, look2)
test.for.zero( look, look3)

set.seed(122)
C<-  rnorm( nrow(x2))
look<- exp(-1*rdist(x1,x2)/4)%*%C
look2<- stationary.cov( x1,x2, theta=4, C=C)
look3<- Exp.cov( x1, x2, theta=4.0, C=C)
test.for.zero( look, look2)
test.for.zero( look, look3)

#### check tranformation of coordinates
V<- matrix( c(2,1,0,4), 2,2)
Vi<- solve( V)
u1<- t(Vi%*% t(x1))
u2<- t(Vi%*% t(x2))

look<- exp(-1*rdist(u1,u2))
look2<- stationary.cov( x1,x2, V= V)
test.for.zero( look, look2)

look<- Wendland(rdist(u1,u2), k=3, dimension=2)
look2<- stationary.cov( x1,x2, V= V, Covariance = "Wendland",
                       k=3, dimension=2)
test.for.zero( look, look2)

### check tapering of covariances
x1<- x[1:5,]
x2<- x[2:6,]
V<- matrix( c(2,1,0,4), 2,2)
Vi<- solve( V)

u1<- x1
u2<- x2

look1a<- exp(-1*rdist(u1,u2))
look1b<-  Wendland(rdist(u1,u2),
                                      k=3, dimension=2, theta= 1)
look1<- look1a*look1b
look2<- stationary.taper.cov( x1,x2, theta=1,
               Taper.args=list( theta=1,k=3, dimension=2), verbose=FALSE)
test.for.zero( look1, as.matrix(look2))


u1<- t(Vi%*% t(x1))
u2<- t(Vi%*% t(x2))


look1a<- exp(-1*rdist(u1,u2))
look1b<-  Wendland(rdist(u1,u2),
                                      k=3, dimension=2, theta= 1.5)
look1<- look1a*look1b
look2<- stationary.taper.cov( x1,x2,V=V,
               Taper.args=list( theta=1.5,k=3, dimension=2), verbose=FALSE)
test.for.zero( look1, as.matrix(look2))


u1<- t(Vi%*% t(x1))
u2<- t(Vi%*% t(x2))


look1a<- Matern(rdist(u1,u2), smoothness=1.5)
look1b<-  Wendland(rdist(u1,u2),
                                      k=3, dimension=2, theta= 1.5)
look1<- look1a*look1b
look2<- stationary.taper.cov( x1,x2,V=V,Covariance=Matern, smoothness=1.5,
               Taper.args=list( theta=1.5,k=3, dimension=2), verbose=FALSE)
test.for.zero( look1, as.matrix(look2))


# some tests of great circle distance


stationary.taper.cov( x[1:3,],x[1:10,] , theta=200, Taper.args= 
       list(k=2,theta=300, dimension=2),
       Dist.args=list( method="greatcircle") )-> temp

# temp is now a tapered 3X10 cross covariance matrix in sparse format. 
# should be identical to
# the direct matrix product

temp2<- Exponential( rdist.earth(x[1:3,],x[1:10,]), range=200) * 
           Wendland(rdist.earth(x[1:3,],x[1:10,]), theta= 300, k=2, dimension=2)

test.for.zero(  as.matrix(temp), temp2, tol=2e-6, tag="taper with great circle")

# example of calling the taper version directly 
# Note that default covariance is exponential and default taper is 
# Wendland (k=2).

stationary.taper.cov( x[1:3,],x[1:10,] , theta=1.5, Taper.args= 
      list(k=2,theta=2.0, dimension=2) )-> temp
# temp is now a tapered 5X10 cross covariance matrix in sparse format. 
# should be identical to
# the direct matrix product

temp2<- Exp.cov( x[1:3,],x[1:10,], theta=1.5) * 
           Wendland(rdist(x[1:3,],x[1:10,]),
                      theta= 2.0, k=2, dimension=2)

test.for.zero(  as.matrix(temp), temp2, tag= "high level test of taper cov")

stationary.taper.cov( x[1:3,],x[1:10,] , range=1.5,
        Taper.args= list(k=2,theta=2.0,
                       dimension=2) )-> temp

test.for.zero(  as.matrix(temp), temp2, tag= "high level test of taper cov")

cat("end tests of V argument in covariances", fill=TRUE)













