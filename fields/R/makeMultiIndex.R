makeMultiIndex<- function(M){
# M are  L integers with product prodM
# Will create a  prodM by L matrix that is all combinations of (1:M[i]) for i =1,2, ...L
# This is organized in the standard array ordering where the first column varies the fastest
#  for M = c( 3,2,4)
# 24 rows     bigIndex =  1,1,1 
#                          2,1,1
#                          3,1,1
#                          1,2,1
#                          2,2,1
#                          3,2,1
#                      etc
#                ending at
#                          2,2,4
#                          3,2,4
                
  L<- length( M)
# qucik return for L =1 and 2  
  if( L ==1){ 
    bigIndex<- cbind(1:M[1] )
  }
  if( L==2){
    bigIndex<-  cbind( rep(1:M[1], M[2]),  
                      rep( 1:M[2], rep(M[1], M[2]) )
    )
  }
  if( L > 2){
    prodM<- prod( M)
    bigIndex<-matrix( NA, prodM,L)
    
    for ( k in 1: L){
      if( k ==1){
        ML<- 1
      }
      else{
        ML<- prod( M[1: (k-1)] )
      }
      
      if( k ==L){
        MU<- 1}
      else{
        MU<-prod( M[ (k+1) :  L  ] )
      }
      gridTmp<- rep( 1:M[k],  rep(ML, M[k]) )
      bigIndex[,k]<- rep( gridTmp, MU)
    }
  }
return( bigIndex )  
  
}
