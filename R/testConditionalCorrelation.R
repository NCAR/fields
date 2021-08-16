
library( fields)
sGrid<- cbind(  rep( 1:5, each=4),rep(1:4,5))

sOff1<- rbind( c(2,2.5),
               c( 3,2.5)
)
sOff2<- rbind( c(3.25,2.5),
               c( 2.75,2.5)
)
pdf("figObsGrid.pdf", width=5, height=4.5)
plot( sGrid, cex=1.25, pch=16, xlab="X", ylab="Y")
xline(1:5 , col="grey", lty=2)
yline( 1:4, col="grey", lty=2)
points( sOff1,  col="red",  pch=16, )
points( sOff2, col="grey", pch=15)
plot( sGrid)
points( sOff, col="red", pch=16)
dev.off()

testCor<- function( sOff, aRange, smoothness){
  Sigma11<- stationary.cov( sGrid,sGrid,
                            Covariance="Matern",
                            aRange= aRange,
                            smoothness=smoothness)
  Sigma21<- stationary.cov( sOff,sGrid, aRange= aRange,
                            Covariance="Matern",
                            smoothness=smoothness)
  Sigma22<- stationary.cov( sOff,sOff, aRange= aRange,
                            Covariance="Matern",
                            smoothness=smoothness)
  covPred<- Sigma22 - Sigma21%*%solve( Sigma11)%*%t(Sigma21)
  cD<- diag(1/sqrt(diag(covPred)))
  corPred<- cD%*%covPred%*%cD
  return(c(corPred[2,1], covPred[1,1] ))
}



M<- 100
aGrid<- seq( .2, 10,,M)
smoothGrid<- seq( .5,1.5,,M)
z<- matrix( NA, M,M)
for( j in 1:M){
  #cat(j, " ")
  for(k  in 1:M){
    z[j,k]<- testCor(sOff1, aGrid[j], smoothGrid[k])
  }
}
out1<- list( x= aGrid, y=smoothGrid, z= abs(z))
#image.plot( out)

quantile( out1$z, c( .9, .95,1.0))
which.max.image(out1)
for( j in 1:M){
  #cat(j, " ")
  for(k  in 1:M){
    z[j,k]<- testCor(sOff2, aGrid[j], smoothGrid[k])
  }
}
out2<- list( x= aGrid, y=smoothGrid, z= abs(z))
#image.plot( out)

quantile( out2$z, c( .9, .95,1.0))
which.max.image(out2)

