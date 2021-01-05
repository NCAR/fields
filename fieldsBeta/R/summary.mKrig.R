summary.mKrig <- function(object, ...) {
  
  outObject<- list()  
  digits<- 4
  summaryObject<- object$summary
  
  if (is.matrix(object$residuals)) {
    n <- nrow(object$residuals)
    nData <- ncol(object$residuals)
  }
  else {
    n <- length(object$residuals)
    nData <- 1
  }
  outObject$nData<-nData
  
  c1 <- "Number of Locations:"
  c2 <- n
  
  if (nData > 1) {
    c1 <- c(c1, "Number of data sets fit:")
    c2 <- c(c2, nData)
    
  }
  c1 <- c(c1, "Degree of polynomial null space ( base model):")
  if(object$m !=0 ){
    c2 <- c(c2, object$m - 1)
  }
  else{
    c2 <- c(c2, NA)
  }
  c1 <- c(c1, "Total number of parameters in base model")
  c2 <- c(c2, object$nt)
  if (object$nZ > 0) {
    c1 <- c(c1, "Number of additional covariates (Z)")
    c2 <- c(c2, object$nZ)
  }
  if (!is.na(object$eff.df)) {
    c1 <- c(c1, " Estimate Eff. degrees of freedom")
    c2 <- c(c2, signif(object$eff.df, digits))
    if (length(object$trA.info) < object$np) {
      c1 <- c(c1, "    Standard Error of Eff. Df ")
      c2 <- c(c2, signif(sd(object$trA.info)/sqrt(length(object$trA.info)), 
                         digits))
    }
  }
  c1 <- c(c1, "Smoothing parameter")
  c2 <- c(c2, signif(summaryObject["lambda"], digits))
  
  if (nData == 1) {
    c1 <- c(c1, "tau (nugget sd)")
    c2 <- c(c2, signif(summaryObject["tau"], digits))
    c1 <- c(c1, "sigma (process sd)")
    c2 <- c(c2, signif(sqrt( summaryObject["sigma2"]), digits))
  }
  
  c1 <- c(c1, "Nonzero entries in covariance")
  c2 <- c(c2, object$nonzero.entries)
  sum <- cbind(c1, c2)
  dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
  ########### save table of information and the call
  outObject$summaryTable<- sum
  outObject$call<- object$call
  outObject$collapseFixedEffect<- object$collapseFixedEffect
  
  ########### information for SE for fixed effects
  if( outObject$collapseFixedEffect | (nData==1)){
    outObject$fixedEffectsCov<- object$fixedEffectsCov
    SE<- sqrt(diag(outObject$fixedEffectsCov))
    beta<-  object$beta[,1]
    print( beta)
    pValue<- pnorm(abs(beta/SE), lower.tail = FALSE)*2
    outObject$fixedEffectsTable<- cbind( signif(beta, digits), 
                                         signif(SE, digits),
                                         signif(pValue, digits)
    )
    if( is.null( object$fixedEffectNames )){
      outObject$fixedEffectNames<- paste0("d",1:(object$nt) )
    }
    else{
      outObject$fixedEffectNames<- object$fixedEffectNames
    }
    dimnames( outObject$fixedEffectsTable) <- list( outObject$fixedEffectNames,
                                                    c("estimate", "SE", "pValue") )
    ########### save covariance information
  } 
  outObject$cov.function<- object$cov.function
  outObject$args<- object$args
  class( outObject)<-"mKrigSummary"
  return( outObject)
}
