timeCEWithNN<- function( m,n,M,np){
  s<- matrix( runif(n*2 ),n,2)
  y<-  rnorm( n)

  #  t1<- system.time(
  # mKrigObject<- mKrig( s,y,cov.function=stationary.cov,
  #                               cov.args= list(aRange=.1,
  #                                       Covariance="Matern",
  #                                         smoothness=1.0),
  #                               lambda=.1
  #                                )
  # )
  
  t1<- system.time(
    mKrigObject<- mKrig( s,y,cov.function= Exp.cov,
                         cov.args= list(aRange=.1),
                         lambda=.1
    )
  )
  
t2<- system.time(
     look<- simLocal.mKrig(mKrigObject, M=M,nx=m, ny=m,
                          gridRefinement = 1, np=np, 
                          giveWarnings=FALSE)
        )
#print( look$timingFull )

  return(
c( m,n,M,
  look$timing, mKrig = t1[3],  draw =t2[3]/M )
  )
}  
  
  
  
  