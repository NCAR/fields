makeBigB<-function(s,m,n, exp_cov){
  #
  # function assumes the grid is 
  # integer locations and 1:m by 1:n
  #  grid and off grid locations need to be transformed to that scale
  # 
  # also assumes the grid extends two cells beyond any off
  # e.g. s  coordinates should be between 
  # 2 and m-3 and 2 and n-3
  #
  M<- nrow( s)
  grid_list <- list(x= 1:m, 
                    y= 1:n)
  class(grid_list)<- "gridList"
  # lower left corner of grid box containing the points
  s0<- trunc( s) - 1
  #print( s0)
  # index when 2D array is unrolled
  s0Index<- s0[,1] + (s0[,2]-1)*m
  # specific to 2nd degree neighborhood
  # 16 points total 
  xshift<- rep( c(0,1,2,3),4)
  yshift<- rep( c(0,1,2,3),each=4 )
  nnXY<- cbind( xshift, yshift)
  #
  # MX16 matrices each row are  row and column indices for 
  # the 16 nn 
  sX<- s0[,1] + matrix( rep( xshift,M),
                        nrow=M, ncol=16, byrow=TRUE)
  if( any( (sX < 1)| (sX>m)) ) {
    stop( "s outside range for X")
  }
  sY<- s0[,2] + matrix( rep( yshift,M),
                        nrow=M, ncol=16, byrow=TRUE)
  if( any( (sY < 1)| (sY>n)) ) {
    stop( "s outside range for Y")
  }
  # indices of all nearest neighbors for unrolled vector.
  # this is an MX16 matrix where indices go from 1 to m*n
  # these work for the unrolled 2D array 
  # 
  sIndex<-  sX + (sY-1)*m
  # differences between nn and the off grid locations
  # for both coordinates
  #print( sX[1,])
  dX<- sX- s[,1]
  dY<- sY- s[,2]
  #print( dX[1,])
  # all pairwise distances between each off grid and 
  # 16 nearest neighbors 
  dAll<- sqrt(dX^2 + dY^2)
  
  # cross covariances 
  Sigma21Star<-  exp_cov(dAll)
  # covariance among nearest neighbors 
  Sigma11 <-  exp_cov( rdist(nnXY, nnXY ) )
 
  Sigma11Inv<- solve( Sigma11)
  #print( Sigma11Inv[1,])
  # each row of B are the weights used to predict off grid point
  B<- Sigma21Star%*%Sigma11Inv
  
  # create spind sparse matrix
  # note need to use unrolled indices to refer to grid points
  ind<- cbind( rep(1:M, each=16), c( t( sIndex)))
  ra<-  c( t( B))
  da<- c( M, m*n )
  spindBigB<-  list(ind=ind, ra=ra, da=da )
  # now convert to the more efficient spam format
  BigB<- spind2spam( spindBigB)
#
#  return( list(BigB= BigB, 
#              sIndex=sIndex,
#              Sigma11Inv= Sigma11Inv,
#              Sigma21Star= Sigma21Star,
#              sX=sX, sY=sY)
#  )

    return(BigB)
  
}
