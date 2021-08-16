

timeOffGrid<- function( m,M, np){
  n<- m
  dx<- 1
  dy<- 1
  sigma2<-2.0 
  
  set.seed( 123)
  # locations random but  avoid edges
  s<- cbind( dx*runif( M,  np, (m-(np+1))),
             dy* runif( M, np, (n-(np+1)))
  )
  
  
  # random uniform is ok as we are just checking agreement
  set.seed( 222)
  y<- matrix( runif(m*n),m,n)
  yUnrolled<- c( y)
  
  #look<- sparseB%*%yUnrolled
  
  ###################################
  # test new function
  ###################################
  t1<- system.time(
    sparseObj0<-  offGridWeights( s, list( x= 1:m, y=1:n),
                                  aRange=10, sigma2=2.0, 
                                  Covariance="Exponential", 
                                  np=np,
                                  giveWarnings=FALSE)
  )
  
  t2<- system.time(
    look5<- sparseObj0$B%*%yUnrolled
  )
  out<-c( m, M, c(t1[3]), c(t2[3]) ) 
  return( out)
}
