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


# test of sreg and related functions

suppressMessages(library(fields))
options(echo=FALSE)
 test.for.zero.flag<- 1

set.seed(123)

# Tps has been tested from scratch using basic linear algebra
# so test sreg against this

x<- rat.diet$t
y<- rat.diet$trt

sreg( x,y, lambda= 10)-> out
Tps( x,y, scale="unscaled", lambda=10*length(y))-> out2

test.for.zero( out$fitted.values, out2$fitted.values, tag="predict at lambda sreg/Tps")


#### GCV test

sreg( x,y, tol=1e-12)-> out
gcv.sreg( out, tol=1e-12)$lambda.est -> look0

test.for.zero( out$lambda.est[1,2], look0[1,2], tol=5e-4)

Tps( x,y)-> out2
gcv.Krig( out2, tol=1e-12)$lambda.est[1,2]-> look2
gcv.sreg( out, tol=1e-12)$lambda.est[1,2] -> look

test.for.zero( look, look2, tol=2.1e-6, tag="GCV sreg/Tps")

#### replications
set.seed( 123)
x<-  rep(rat.diet$t,3)
y<- rep( rat.diet$trt,3) + rnorm(39*3)*5

sreg( x,y)-> out
gcv.sreg( out, tol=1e-12)$lambda.est -> look

Tps( x,y, scale="unscaled")-> out2
gcv.Krig( out2, tol=1e-12)$lambda.est-> look2
look2[,1]<- look2[,1]/length( out$xM)

test.for.zero( look[1:3,3], look2[1:3,3],  
              tag="GCV sreg/Tps reps case",tol=1e-06)

test.for.zero( look[2,3], look2[2,3], tol=1e-6, 
              tag="GCV sreg/Tps reps case")


cat( "All done with sreg tests", fill=TRUE)
options(echo=TRUE)









