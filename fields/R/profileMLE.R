profileMLE<- function (obj, parName, parGrid=NULL, gridN=15,
                       cov.params.start=NULL, GCV=FALSE, REML=FALSE,
                       verbose=FALSE ){
  if( class(obj)[1]!= "spatialProcess"){
    stop("only implemented tfor spatialProcess objects")
  }
  if( obj$CASE==0){
    stop("object needs to be have MLE parameters")
  }
  if(is.null(parGrid) ){
  parGrid<- data.frame(obj$summary[parName] *
                         seq( .5,2,length.out=gridN)
                       )
  names(parGrid)<- parName
  }
  
  # remove profile parameter as having a starting value
  
  if( is.null(cov.params.start)){
  cov.params.startTmp<-  as.list(obj$MLEInfo$pars.MLE)
  #cov.params.startTmp<- obj$cov.params.start
  cov.params.startTmp[parName]<- NULL
  }
  else{
    cov.params.startTmp<-cov.params.start
  }
  
# just evaluate other parameters at MLEs and don't optimze. 
# E.g. just use  MLE lambda hat for all choice in in parGrid. 
# if( fast){
#    cov.params.startTmp<-NULL
#    
#  }
#  cov.argsTemp<- obj$cov.argsFull
#  cov.argsTemp[parName]<-  NULL
  
  profileInfo<-  mKrigMLEGrid(obj$x, obj$y,  
               weights = obj$weights,
               Z = obj$Z, 
               mKrig.args = obj$mKrig.args,
               cov.function = obj$cov.function, 
               cov.args  = obj$cov.args,
               par.grid = parGrid, 
               na.rm = TRUE,
               verbose = verbose,
               REML = REML,
               GCV  =  GCV,
               cov.params.start = cov.params.startTmp)
  return( profileInfo)
}