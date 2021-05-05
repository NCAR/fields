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
  
  c1 <- c(c1, "tau  Nugget stan. dev:")
  c2 <- c(c2, signif(object$MLESummary["tau"], digits))
    
  c1 <- c(c1, "sigma Process variance: ")
  c2 <- c(c2, signif( sqrt(object$MLESummary["sigma2"]), digits))
  
  c1 <- c(c1, "lambda   tau^2/sigma^2: ")
  c2 <- c(c2, signif(object$MLESummary["lambda"], digits))
  
  c1 <- c(c1, "aRange parameter (in units of distance): ")
  c2 <- c(c2, signif(object$MLESummary["aRange"], digits))
  
  if (!is.na(object$eff.df)) {
    c1 <- c(c1, "Approx.  degrees of freedom for curve")
    c2 <- c(c2, signif(object$eff.df, digits))
    if (length(object$trA.info) < object$np) {
      c1 <- c(c1, "   Standard Error of df estimate: ")
      c2 <- c(c2, signif(sd(object$trA.info)/sqrt(length(object$trA.info)), 
                         digits))
    }
  }
  
  c1<- c(c1, "log Likelihood: " )
  c2<- c( c2, object$summary["lnProfileLike.FULL"])
  c1<- c(c1, "log Likelihood REML: " )
  c2<- c( c2, object$summary["lnProfileREML.FULL"])
  
  summaryStuff<-  cbind(c1, c2)
  dimnames(summaryStuff) <- list(rep("",
                               dim(summaryStuff)[1]), 
                               rep("", dim(summaryStuff)[2]))
###########  
  outObject$summaryTable<- summaryStuff
  outObject$collapseFixedEffect<- object$collapseFixedEffect
  outObject$CITable<- object$CITable
###########
  if( !is.null( object$MLEInfo)){
  outObject$MLEpars<-  names( object$MLEInfo$pars.MLE) 
  outObject$MLESummary<- object$summary
  }
  else{
    outObject$MLEpars <- NA
    outObject$MLESummary<- object$summary
  }
########### information for SE for fixed effects
  if(!is.null(object$beta) ){
  if( outObject$collapseFixedEffect | (nData==1) ){
    outObject$fixedEffectsCov<- object$fixedEffectsCov
    SE<- sqrt(diag(outObject$fixedEffectsCov))
    beta<-  object$beta[,1]
    pValue<- pnorm(abs(beta/SE), lower.tail = FALSE)*2
    outObject$fixedEffectsTable<- cbind( signif(beta, digits), 
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
  }
  else{
    outObject$fixedEffectsTable<- NA
  }
#####################
  outObject$nData <- nData
  outObject$call<- object$call
  outObject$cov.function<- object$cov.function
  outObject$args<- object$args
  outObject$nonzero.entries<- object$nonzero.entries
  outObject$MLEInfo<- object$MLEInfo
  
  class( outObject)<-"spatialProcessSummary"
    
  return( outObject)

}
