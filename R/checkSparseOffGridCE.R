# Local Kriging - sparse matrix implementation (small example)
setwd("~/Dropbox/Home/Repositories/fields")
rm(list=ls())

library(fields)
library(LatticeKrig)
source( "makeBigB.R")
source("offGridWeights.R")

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


#s0<- trunc( s) - 1
colTab<- tim.colors(M)
plot( s, xlim=c(1,m), ylim= c(1,n),
       pch=16, col="magenta", cex=1.5)
xline( 1:m, col="grey")
yline( 1:n, col="grey")

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
    points( sTemp, pch=16, col=colTab[k])
  }
  Sigma11<- sigma2*exp_cov(rdist(sTemp, sTemp))
  Sigma11Inv<- solve(Sigma11)
  Sigma21<- sigma2*exp_cov(rdist(rbind(s[k,]), sTemp))
  Btest<- Sigma21%*%Sigma11Inv
  result<- Sigma21%*%Sigma11Inv%*%yTemp
  look2[k]<- result
  predVar[k]<- sigma2 - diag(Sigma21%*%Sigma11Inv%*%t(Sigma21 ) )
}

points( s, pch=16, col="magenta")
points( s,  col="black")

if( np==2){
# only works for np=2
sparseB<- makeBigB( s, m,n, exp_cov)
look<- sparseB%*%yUnrolled
test.for.zero( look, look2)
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
fakeObj<- list( args = list( Covariance= "Exponential" ),               ),
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

















