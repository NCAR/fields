profileCI<- function( obj, parName, 
                      CIlevel=.95, 
                      REML=FALSE){
  y<- obj$summary[,ifelse(REML,"lnProfileREML.FULL",
                           "lnProfileLike.FULL")]
  x<- obj$summary[,parName]
  logX<- log10( x)
  lGrid<- (seq( min(logX), max(logX),length.out=400))
  yGrid<- splint( logX, y, lGrid)
  yMax<- max(yGrid)
  ChiLevel<- qchisq(.95,1)
  ind<-  yMax - yGrid <= ChiLevel 
  CI<- range( lGrid[ind] )
  CI <- 10**CI
  return( CI)
}