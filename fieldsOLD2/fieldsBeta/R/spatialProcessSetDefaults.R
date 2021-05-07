spatialProcessSetDefaults<- function( cov.function,
                                      cov.args,
                                      cov.params.start,
                                      mKrig.args,
                                      extraArgs=NULL,
                                      profileLambda,
                                      profileARange,
                                      gridARange,
                                      gridLambda,
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
  
  # overwrite the default choices is some are passed as ...
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
 
 
  
  #CASE 0 is to evaluate at fixed lambda and aRange
  if( !is.null( cov.args$lambda) & !is.null( cov.args$aRange)){
    CASE<- 0
  }
  
  if( is.null( cov.args$lambda) & is.null(cov.params.start$lambda) ){
    cov.params.start$lambda <- .5
  }
  
  #CASE 1 is to find MLE for  lambda but aRange is fixed
  # some lambda info
  if( !is.null(cov.params.start$lambda)& !is.null( cov.args$aRange) ){
    
    CASE<- 1
  }
  
  # CASE 2 is to find MLEs for both lambda and aRange
  # a grid search is done in place of a single starting value for aRange
  
  if( is.null( cov.args$aRange) & !is.null(cov.params.start$aRange) ){
    CASE<- 2
  }
  
  # CASE 3 is to find MLEs for both lambda and aRange
  # a grid search is done in place of a single starting value for aRange
  # the nuclear option ...
  noARange<- is.null( cov.args$aRange) & is.null(cov.params.start$aRange) 
  if( noARange | profileARange ){
    cov.params.start$aRange <- NA
    CASE<- 3
  }
  
  # CASE 5 is to find MLEs for both lambda and aRange
  # a grid search is done in place of a single starting value for aRange
  # and profiling is done on lambda using grid serach on aRange. 
  if( CASE ==3 & (profileLambda | !is.null(gridLambda) )  ){
    CASE<- 5
  }
  
 # CASE 4 has fixed aRange but profiling/ grid search on lambda 
  if( CASE==1 & (profileLambda | !is.null(gridLambda) )  ){
    CASE<- 4
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
  if( is.null(cov.params.start)){
    cov.params.start<- list( lambda =  .5) 
  }
  
  if( is.null(cov.params.start$lambda) ) {
    cov.params.start$lambda <-  .5 
  }
  
  
  out<- 
    list(  
        cov.function = cov.function,
            cov.args = cov.args,
    cov.params.start = cov.params.start,
          mKrig.args = mKrig.args, 
                CASE = CASE
        )
  
 
 
  return(
         out
          )
}