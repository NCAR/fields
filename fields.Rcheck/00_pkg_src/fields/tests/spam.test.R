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


# test of rdist.near


suppressMessages(library(fields))
options(echo=FALSE)
test.for.zero.flag<- 1

set.seed(123)
x1<- matrix( runif(2*20), ncol=2)
x2<- matrix( runif(2*10), ncol=2)

fields.rdist.near( x1,x2, delta=.75)-> look
temp<- matrix( NA, nrow(x1),nrow(x2))
temp[ look$ind] <- look$ra
temp2<- rdist( x1, x2)
temp2[ temp2> .75] <- NA
temp[ is.na( temp)]<- 0
temp2[ is.na( temp2)]<- 0

test.for.zero( temp, temp2)


# test of constructing covariance matrix
# and also versions of Wendland function
# default taper is wendland k=2.
DD<- rdist( x1,x2)
temp<- Wendland2.2(DD, aRange=.8)
temp2<- Wendland( DD, aRange=.8, dimension=2, k=2)

test.for.zero( temp, temp2)




stationary.taper.cov( x1,x2, Taper="Wendland2.2", 
           Taper.args= list( aRange=.8), spam.format=FALSE )-> look
temp0<- look

stationary.taper.cov( x1,x2, Taper="Wendland2.2",
           Taper.args= list( aRange=.8), spam.format=TRUE )-> look
temp1<-  spam2full( look)

test.for.zero( temp1, temp0)

stationary.taper.cov( x1,x2, Taper="Wendland",
           Taper.args= list( aRange=.8, k=2, dimension=2),
                     spam.format=TRUE )-> look
temp1b<-  spam2full( look)


temp2<-  Wendland2.2(DD, aRange=.8) * Exponential(DD)
temp3<-  wendland.cov(x1,x2, k=2, aRange=.8) * Exponential(DD)
temp4<-  Wendland(DD, k=2, dimension=2, aRange=.8)* Exponential(DD)


test.for.zero( temp1, temp0, rel=FALSE)
test.for.zero( temp1b, temp0, rel=FALSE)
test.for.zero( temp2, temp0, rel=FALSE)

test.for.zero( temp2, temp3,rel=FALSE)
test.for.zero( temp2, temp4,rel=FALSE)



set.seed( 256)
rv<- runif( nrow(x2))

# test of multiply 
stationary.taper.cov( x1, x2,  C= rv)-> look
temp2<-stationary.taper.cov( x1,x2)

spam2full(temp2)%*%(rv)-> look2
test.for.zero( look, look2)

#

set.seed( 123)
temp<- matrix( 1:48, ncol=6, nrow=8)
temp[ sample( 1:48, 20)] <- 0

as.spam( temp)-> temp2
test.for.zero( spam2full(temp2), temp )

spam2spind( temp2)-> temp3

test.for.zero( spind2full( temp3), temp)

test.for.zero( spind2spam( temp3),temp2)

# test that ordering works
MM<- nrow( temp3$ind)
ix<-  sample( 1:MM,MM)
# shuffle temp3
temp3$ind<- temp3$ind[ix,]
temp3$ra<- temp3$ra[ix]

test.for.zero( spind2spam( temp3),temp2)



# temp<- temp[1:4, 1:5]  for help file
#

set.seed( 234)

CC<- matrix( rnorm( 64), 8,8)
A<- ( CC)%*% t(CC)
as.spam( A)-> As

test.for.zero( solve( As), solve( A))

set.seed( 233)
CC<- diag( 1, 8)
CC[4,1:8] <- rnorm(8)
CC[7,1:8] <- rnorm(8)
A<- ( CC)%*% t(CC)
as.spam( A)-> As

test.for.zero( solve( As), solve( A))


data( ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]


Krig(x,y, cov.function = "stationary.taper.cov", aRange=1.5,
      give.warnings=FALSE, 
      cov.args= list( spam.format=FALSE, 
           Taper.args= list( dimension=2, aRange=2.0,k=3) )    ) -> out1

Krig(x,y, cov.function = "stationary.taper.cov", lambda=2.0, aRange=1.5,
      cov.args= list( spam.format=TRUE,
        Taper.args= list( aRange=2.0,k=3, dimension=2)  )
           ) -> out2

temp1<- predict( out1,lambda=2.0)
temp2<- predict( out2)
test.for.zero( temp1, temp2)

cat( "All done with SPAM tests", fill=TRUE)
options(echo=TRUE)
