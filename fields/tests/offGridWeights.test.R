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

# Local Kriging - sparse matrix implementation (small example)

#source( "makeBigB.R")

# simple covariance function for implementation
exp_cov <- function(dist){
  sigma2<-1 
  covariance <- sigma2* exp(-dist / 10) # 10 is arbitrary 
  return(covariance)
}




# -----------------------------
# Define grid and observations
# -----------------------------
 
m<- 40
n<- 45
nx<- m
ny<- n
M<- 10
dx<- 1
dy<- 1
sigma2<-2.0 
np<-3 

set.seed( 123)
# locations random but  avoid edges
s<- cbind( dx*runif( M,  np, (m-(np+1))),
           dy* runif( M, np, (n-(np+1)))
)


# random uniform is ok as we just checking agreement
set.seed( 222)
y<- matrix( runif(m*n),m,n)
yUnrolled<- c( y)

#look<- sparseB%*%yUnrolled



look2<- predVar<- rep( NA, M)
theShift<- 0:(2*np-1) - (np - 1) 
for(  k in 1:M){
  #k<- 3
  yTemp<- NULL
  sTemp<- NULL
  i0<- trunc(s[k,1])
  j0<- trunc(s[k,2])
  for( j in theShift + j0){
     
     for( i in theShift + i0){
       #cat( i,j, fill=TRUE)
    sTemp<- rbind( sTemp, c(i,j))
    yTemp<-  c(yTemp,y[i,j])
     }
   
  }
  Sigma11<- sigma2*exp_cov(rdist(sTemp, sTemp))
  Sigma11Inv<- solve(Sigma11)
  Sigma21<- sigma2*exp_cov(rdist(rbind(s[k,]), sTemp))
  Btest<- Sigma21%*%Sigma11Inv
  result<- Sigma21%*%Sigma11Inv%*%yTemp
  look2[k]<- result
  predVar[k]<- sigma2 - diag(Sigma21%*%Sigma11Inv%*%t(Sigma21 ) )
}



###################################
# test new function
###################################
sparseObj0<-  offGridWeights( s, list( x= 1:m, y=1:n),
                              aRange=10, sigma2=sigma2, 
                              Covariance="Exponential", 
                              np=np)
look5<- sparseObj0$B%*%yUnrolled

test.for.zero( look2, look5 )

test.for.zero(predVar, sparseObj0$predictionVariance )


mKrigObj<- mKrig( s, rnorm( nrow(s)),
                            sigma2=sigma2, tau=0,
                            aRange=10)

sparseObj<-  offGridWeights( s, list( x= 1:m, y=1:n),
                mKrigObject = mKrigObj, np=np
                )

look3<- sparseObj$B%*%yUnrolled

test.for.zero( look2, look3 )

test.for.zero(predVar, sparseObj$predictionVariance )

# cheating on mKrig object
fakeObj<- list( args = list( Covariance= "Exponential" ), 
                summary= c(aRange=10*2.5, sigma2=sigma2)
                )
sparseObj1<-  offGridWeights( s*2.5,
                              list( x = (1:m)*2.5,
                                    y = (1:n)*2.5 ),
                              mKrigObject = fakeObj,
                                       np = np
)

look4<- sparseObj1$B%*%yUnrolled

test.for.zero( look2, look4 )

test.for.zero( sparseObj$predictionVariance, 
              predVar )

















