spatialProcessSetDefaults<- function( cov.function,
                                      cov.args,
                                      cov.params.start,
                                      mKrig.args,
                                      extraArgs=NULL,
                                      parGrid,
                                      gridN,
                                      doGridSearch,
                                      verbose=FALSE)
{
  
  ## convenient defaults for GP fitting.
  ## and also sort what starting parameter values are provided
  #
  # aRange and lambda  are  handled specially because are almost always 
  # estimated and this will simplify the call in  this top level function 
  #
  if( is.null( cov.function)){
    cov.function <- 'stationary.cov'
    if( is.null(cov.args) ){
      cov.args<- list()
    }
    if( is.null(cov.args$Covariance )){
      cov.args$Covariance<- "Matern"
      if( is.null(cov.args$smoothness )){
        cov.args$smoothness<- 1.0
      }
    }
  } 
  
  # overwrite the default choices if some are passed as ...
  #  (some R arcania!)
  
  if( !is.null( extraArgs)){
    if(!is.null(cov.args)){
      ind<- match( names(cov.args), names(extraArgs) ) 
      cov.args <- c( cov.args[is.na(ind)], (extraArgs) )
    }
    else{
      cov.args <- list(extraArgs)
    }
  }
 
  if( verbose){
    cat("update and passed cov.args", fill=TRUE)
    print( cov.args)
  }
  
  # CASE 0 is to evaluate at fixed lambda and aRange
  # and there are no other parameters to optimize over.
  
  if( !is.null( cov.args$lambda) & 
      !is.null( cov.args$aRange) &
       is.null( cov.params.start) 
       ){
    CASE<- 0
  }
  
  #CASE 1 is to find MLEs using starting values provided a grid has not been 
  # supplied for an initial grid search.
  
  if( !is.null(cov.params.start) & is.null(parGrid) ){
    CASE<- 1
  }
 
  if( !is.null(parGrid) ){
    CASE<- 2
  }
  
  noLambda<- is.null( cov.args$lambda) & is.null(cov.params.start$lambda)
  noARange<- is.null( cov.args$aRange) & is.null(cov.params.start$aRange)
 
  if( noLambda |  noARange ){
  print( "create parGrid")
  }
  
# linear fixed model if not specified. 
  if( is.null(mKrig.args)){
    mKrig.args<- list( m=2)
  }
  
# don't find eff df for optimization  
  if( is.null(mKrig.args$find.trA) ){
    if( (CASE >=3)){
    mKrig.args<- c( mKrig.args, list(find.trA = FALSE))
  }
  else{
    mKrig.args<- c( mKrig.args, list(find.trA = TRUE))
  }
  }
  
  #
  # tuck in starting value for lambda if missing
  # 
 
  
  out<- 
    list(  
        cov.function = cov.function,
            cov.args = cov.args,
          mKrig.args = mKrig.args, 
                CASE = CASE
        )
  
 
 
  return(
         out
          )
}