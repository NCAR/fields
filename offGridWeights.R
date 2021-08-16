offGridWeights<-function(s, gridList, np=2,
                         mKrigObject=NULL, 
                         Covariance=NULL, covArgs=NULL,
                         aRange=NULL, sigma2=NULL
                   ){
  #
  # function assumes the grid is 
  # integer locations and 1:m by 1:n
  # grid and off grid locations need to be transformed to that scale
  # 
  # also assumes the grid extends two cells beyond any off
  # e.g. s  coordinates should be between 
  # 2 and m-3 and 2 and n-3
  #
  # If mKrigObject (result of fitting model) is given 
  # extract all the covariance information from it. 
  if( !is.null( mKrigObject)){
    sigma2<- mKrigObject$summary["sigma2"]
    aRange<- mKrigObject$summary["aRange"]
    Covariance<- mKrigObject$args$Covariance
    if( is.null(Covariance)){
      Covariance<- "Exponential"
    }
    covArgs<-mKrigObject$args 
  # some R arcania -- strip out all arguments used by say stationary.cov
  # but not used by the Covariance function 
  # Typically for the Matern family all that is left is smoothness
    if( !is.null( covArgs) ){
      argNames<- names( as.list( get(Covariance)))
      argNames<- argNames[ -length(argNames)]
      ind<- match(  names(covArgs), argNames)
      covArgs[is.na( ind)] <- NULL
    }
  }
  
  m<- length( gridList$x)
  n<- length( gridList$y)
  
  dx<- gridList$x[2]- gridList$x[1]
  dy<- gridList$y[2]- gridList$y[1]
  
  M<- nrow( s)
  # lower left corner of grid box containing the points
  s0<-  cbind( 
               trunc( (s[,1]- gridList$x[1] )/dx) + 1 ,
               trunc( (s[,2]- gridList$y[1] )/dy) + 1
               ) 
  
  # index when 2D array is unrolled
  s0Index<- s0[,1] + (s0[,2]-1)*m
  # np=2
  # specific to 2nd degree neighborhood
  #   (2*np)^2 = 16  points total 
  #xshift<- rep( c(0,1,2,3), 4)
  #yshift<- rep( c(0,1,2,3), each=4 )
  
  theShift<- (0:(2*np-1)) - (np-1)
  xshift<- rep( theShift, 2*np)
  yshift<- rep( theShift, each=2*np )
  
  nnXY<- cbind( xshift, yshift)
  nnXYCoords<- cbind( xshift*dx, yshift*dy)
  
  #
  #  sX and sY are M by (2*np)^2 matrices  where each row is
  #  the unrolled row and column indices for  np=2
  #  yield 16 nearest neighbors 
  
  sX<- s0[,1] + matrix( rep( xshift,M),
                        nrow=M, ncol=(2*np)^2, byrow=TRUE)
  if( any( (sX < 1)| (sX>m)) ) {
    stop( "s outside range for grid")
  }
  sY<- s0[,2] + matrix( rep( yshift,M),
                        nrow=M, ncol=(2*np)^2, byrow=TRUE)
  if( any( (sY < 1)| (sY>n)) ) {
    stop( "s outside range for grid")
  }
  # indices of all nearest neighbors for unrolled vector.
  # this is an M by (2*np)^2 matrix where indices go from 1 to m*n
  # these work for the unrolled 2D array 
  # 
  sIndex<-  sX + (sY-1)*m
  # differences between nn and the off grid locations
  # for both coordinates
  # convert from integer grid to actual units. 
  differenceX<- (sX-1)*dx + gridList$x[1] - s[,1]
  differenceY<- (sY-1)*dy + gridList$y[1] - s[,2]
  # all pairwise distances between each off grid and 
  # (2*np)^2  ( np=2 has 16) nearest neighbors 
  dAll<- sqrt(differenceX^2 + differenceY^2)
  # pairwise distance among nearest neighbors. 
  dNN<- rdist(nnXYCoords, nnXYCoords )
  # cross covariances
  Sigma21Star<- sigma2* do.call(Covariance,
                                c(list(d = dAll/aRange), 
                                         covArgs)) 
  # covariance among nearest neighbors 
  Sigma11 <-  sigma2* do.call(Covariance,
                              c(list(d = dNN/aRange), 
                                covArgs))
  Sigma11Inv <- solve( Sigma11)
  # each row of B are the weights used to predict off grid point
  B <- Sigma21Star%*%Sigma11Inv
  # prediction variances  
  # use cholesky for more stable numerics
  cholSigma11Inv<- chol(Sigma11Inv)
  #  sigma2 - diag(Sigma21Star%*%Sigma11Inv%*%t(Sigma21Star) )
  w <- Sigma21Star%*%t(cholSigma11Inv)
  predictionVariance <-  sigma2 - rowSums(w^2)
  # create spind sparse matrix
  # note need to use unrolled indices to refer to grid points
  ind<- cbind( rep(1:M, each= (2*np)^2 ), c( t( sIndex)))
  ra<-  c( t( B))
  da<- c( M, m*n )
  spindBigB<-  list(ind=ind, ra=ra, da=da )
  # now convert to the more efficient spam format
  BigB<- spind2spam( spindBigB)
    return(
      list( B= BigB, predictionVariance = predictionVariance)
          )
  
}
