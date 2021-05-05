confidenceIntervalMLE<- function( obj, CILevel, verbose=FALSE){
   if( is.na( CILevel)){
     cat("CI not found")
     return(NULL)
   }
  MLEInfo<- obj$MLEInfo
  transformedMLEPars<- MLEInfo$optimResults$par
  sampleFisherInfo<- solve(-1*MLEInfo$optimResults$hessian)
  if( verbose){
    print( sampleFisherInfo)
  }
  SE<-  sqrt( diag(sampleFisherInfo ))
  CITmp1<- transformedMLEPars - qnorm(CILevel)*SE
  CITmp2<- transformedMLEPars + qnorm(CILevel)*SE
  tableCI<- 
    cbind(
      MLEInfo$parTransform( CITmp1, inv=TRUE),
      MLEInfo$parTransform( CITmp2, inv=TRUE)
    )
  percentCI<- paste0(100* CILevel,"%")
   colnames(tableCI ) <- c(
    paste0("lower",percentCI), 
    paste0("upper",percentCI)
  )
  return( tableCI)
}