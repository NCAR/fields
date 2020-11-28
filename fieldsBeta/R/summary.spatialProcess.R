# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2018
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
summary.spatialProcess <- function(object, ...) {
# output list
  outObject<- list()
  digits<- 4
  if (is.matrix(object$residuals)) {
    n <- nrow(object$residuals)
    nData <- ncol(object$residuals)
  }
  else {
    n <- length(object$residuals)
    nData <- 1
  }
  
  c1 <- "Number of Observations:"
  c2 <- n
  
  if (nData > 1) {
    c1 <- c(c1, "Number of data sets fit:")
    c2 <- c(c2, nData)
  }
  
  c1 <- c(c1, "Degree of polynomial in fixed part: ")
  
  
  if(object$m !=0 ){
    c2 <- c(c2, object$m - 1)
  }
  else{
    c2 <- c(c2, NA)
  }
  c1 <- c(c1, "Total number of parameters in fixed part: ")
  c2 <- c(c2, object$nt)
  if (object$nZ > 0) {
    c1 <- c(c1, "Number of additional covariates (Z)")
    c2 <- c(c2, object$nZ)
  }
  
  c1 <- c(c1, "MLE nugget variance ( tau^2)")
  c2 <- c(c2, signif(object$tau.MLE.FULL^2, digits))
    
  c1 <- c(c1, "MLE process variance (sigma)")
  c2 <- c(c2, signif(object$sigma.MLE.FULL, digits))
  
  c1 <- c(c1, "Value for lambda = tau^2/sigma")
  c2 <- c(c2, signif(object$lambdaModel, digits))
  
  
  c1 <- c(c1, "MLE range parameter (theta, units of distance): ")
  c2 <- c(c2, signif(object$theta.MLE, digits))
  
  c1 <- c(c1, paste0( "Approx ", object$confidenceLevel,  "% CI for theta:  ") ) 
  c2<-  c(c2, paste( "[",signif( object$theta.CI[1], digits), ",",
                      signif( object$theta.CI[2], digits), "]"  ) 
           )
  
  if (!is.na(object$eff.df)) {
    c1 <- c(c1, "Approx.  degrees of freedom for curve")
    c2 <- c(c2, signif(object$eff.df, digits))
    if (length(object$trA.info) < object$np) {
      c1 <- c(c1, "   Standard Error of df estimate: ")
      c2 <- c(c2, signif(sd(object$trA.info)/sqrt(length(object$trA.info)), 
                         digits))
    }
  }
  
  c1 <- c(c1, "Nonzero entries in covariance")
  c2 <- c(c2, object$nonzero.entries)
  
  c1<- c(c1, "log Likelihood: " )
  c2<- c( c2, object$lnProfileLike.FULL)
  c1<- c(c1, "log Likelihood REML: " )
  c2<- c( c2, object$lnProfileREML.FULL)
  
  summaryStuff<-  cbind(c1, c2)
  dimnames(summaryStuff) <- list(rep("",
                               dim(summaryStuff)[1]), 
                               rep("", dim(summaryStuff)[2]))
###########  
  outObject$summaryTable<- summaryStuff
  outObject$collapseFixedEffect<- object$collapseFixedEffect
###########
  outObject$MLEpars<-  names( object$MLEInfo$pars.MLE) 
  outObject$MLESummary<- object$summary
########### information for SE for fixed effects
  if( outObject$collapseFixedEffect | (nData==1) ){
    outObject$fixedEffectsCov<- object$fixedEffectsCov
    SE<- sqrt(diag(outObject$fixedEffectsCov))
    d.coef<-  object$d[,1]
    pValue<- pnorm(abs(d.coef/SE), lower.tail = FALSE)*2
    outObject$fixedEffectsTable<- cbind( signif(d.coef, digits), 
                                      signif(SE, digits),
                                      signif(pValue, digits)
                                    )
    if( is.null( object$fixedEffectNames ) ){
      outObject$fixedEffectNames<- paste0("d",1:(object$nt) )
    }
    else{
      outObject$fixedEffectNames<- object$fixedEffectNames
    }
    dimnames( outObject$fixedEffectsTable) <- list( outObject$fixedEffectNames,
                                           c("estimate", "SE", "pValue") )
  }
#####################
  outObject$nData <- nData
  outObject$call<- object$call
  outObject$cov.function<- object$cov.function
  outObject$args<- object$args
 
    class( outObject)<-"spatialProcessSummary"
    
  return( outObject)

}
