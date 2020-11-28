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
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    

mKrig.trace <- function(object, iseed, NtrA) {
    set.seed(iseed)
    # if more MonteCarlo samples > number of data points just
    # find A exactly using  np  calls to predict.
    np<- object$np
    if (NtrA >= object$np) {
        Ey <- diag(1, np)
        NtrA <- np
        hold <- diag(predict.mKrig(object, ynew = Ey, collapseFixedEffect=FALSE))
        trA.info<- NA
        trA.est <- sum(hold)
    }
    else {
        # if fewer tests then use random trace method
        # find fitted.values  for iid N(0,1) 'data' to calculate the
        # the Monte Carlo estimate of tr A(lambda)
        # basically repeat the steps above but take some
        # short cuts because we only need fitted.values
        # create random normal 'data'
        Ey <- matrix(rnorm(np * NtrA), nrow = np, 
            ncol = NtrA)
        trA.info <- colSums(Ey * (predict.mKrig(object, ynew = Ey,
                                    collapseFixedEffect=FALSE)))
        trA.est <- mean(trA.info)
    }
    if (NtrA < np) {
     MSE<-(sum(object$residuals^2)/np) 
     GCV <-       MSE/(1 - trA.est /np)^2
     GCV.info <- MSE/( 1 - trA.info/np)^2
    }
    else{
    	GCV<- NA
    	GCV.info <- NA
    }	
    return(
    list(trA.info = trA.info, eff.df = trA.est,
             GCV= GCV, GCV.info=GCV.info)
    )
}

mKrig.coef <- function(object, y, collapseFixedEffect=TRUE) {
    # given new data y and the matrix decompositions in the
    # mKrig object find coefficients d and c.
    # d are the coefficients for the fixed part
    # in this case hard coded for a low order polynomial
    # c are coefficients for the basis functions derived from the
    # covariance function.
    #
    # see mKrig itself for more comments on the linear algebra
    #
    # Note that all these expressions make sense if y is a matrix
    # of several data sets and one is solving for the coefficients
    # of all of these at once. In this case d.coef and c.coef are matrices
    #
    # generalized least squares for d
    if( any(is.na(y))){
    	stop("mKrig can not omit missing values in observation vecotor")
    }
   if( object$nt>0){
    d.coef <- as.matrix(qr.coef(object$qr.VT, forwardsolve(object$Mc, 
        transpose = TRUE, y, upper.tri = TRUE)))
    d.coefMeans<- rowMeans( d.coef)
    if( collapseFixedEffect){ 
      d.coef<- matrix( d.coefMeans, ncol=ncol(d.coef), nrow= nrow( d.coef))
    }
    #  residuals from subtracting off fixed part
    #  of model as m-1 order polynomial
    resid <- y - object$Tmatrix %*% d.coef
   }
  else{
    d.coef<- NULL
    resid <- y
  }
    # and now find c.
    c.coef <- forwardsolve(object$Mc, transpose = TRUE, resid, 
        upper.tri = TRUE)
    c.coef <- as.matrix(backsolve(object$Mc, c.coef))
    out <- list(d = (d.coef), c = (c.coef))
    return(out)
}

summary.mKrig <- function(object, ...) {
    
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
    c1 <- c(c1, " Eff. degrees of freedom")
    c2 <- c(c2, signif(object$eff.df, digits))
    if (length(object$trA.info) < object$np) {
      c1 <- c(c1, "   Standard Error of estimate: ")
      c2 <- c(c2, signif(sd(object$trA.info)/sqrt(length(object$trA.info)), 
                         digits))
    }
  }
  c1 <- c(c1, "Smoothing parameter")
  c2 <- c(c2, signif(object$lambda.fixed, digits))
  
  if (nData == 1) {
    c1 <- c(c1, "MLE sigma ")
    c2 <- c(c2, signif(object$shat.MLE, digits))
    c1 <- c(c1, "MLE rho")
    c2 <- c(c2, signif(object$rho.MLE, digits))
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
   d.coef<-  object$d[,1]
   print( d.coef)
   pValue<- pnorm(abs(d.coef/SE), lower.tail = FALSE)*2
   outObject$fixedEffectsTable<- cbind( signif(d.coef, digits), 
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

print.mKrig<- function (x, digits = 4,...){
  object<- summary.mKrig( x,...)
  print.mKrigSummary(object, digits = digits)
}

print.mKrigSummary <- function (x, digits = 4, ...) 
  {
    cat("Call:\n")
    dput(x$call)
    print(x$summaryTable, quote = FALSE)
    
# fixed effects are reported differently when fields are replicated.    
    nData<- x$nData
    cat(" ", fill = TRUE)
    if( nData == 1 | x$collapseFixedEffect ){
      cat(" ", fill = TRUE)
      cat("Summary of fixed effects", fill = TRUE)
      print( x$fixedEffectsTable)
      }
      else {
        cat("Estimated fixed effects found separately for each replicate field", 
            fill = TRUE)
      }
    cat(" ", fill = TRUE)
    cat("Covariance Model:", x$cov.function, fill = TRUE)
    if (x$cov.function == "stationary.cov") {
      cat("   Covariance function:  ", ifelse(is.null(x$args$Covariance), 
                                              "Exponential", x$args$Covariance), fill = TRUE)
    }
    if (!is.null(x$args)) {
      cat("   Non-default covariance arguments and their values ", 
          fill = TRUE)
      nlist <- as.character(names(x$args))
      NL <- length(nlist)
      for (k in 1:NL) {
        cat("   Argument:", nlist[k], " ")
        if (object.size(x$args[[k]]) <= 1024) {
          cat("has the value(s): ", fill = TRUE)
          print(x$args[[k]])
        }
        else {
          cat("too large to print value, size > 1K ...", 
              fill = TRUE)
        }
      }
    }
    invisible(x)
  }

predict.mKrig <- function(object, xnew = NULL, ynew = NULL, grid.list=NULL,
    derivative = 0, Z = NULL, drop.Z = FALSE, just.fixed = FALSE,
    collapseFixedEffect = object$collapseFixedEffect, 
    ...) {
    # the main reason to pass new args to the covariance is to increase
    # the temp space size for sparse multiplications
    # other optional arguments that typically describe the covariance function 
    # from mKrig are passed along in the list object$args
    cov.args <- list(...)
    # predict at observation locations by default
    if( !is.null(grid.list)){
        xnew<- make.surface.grid(grid.list)
      }
    if (is.null(xnew)) {
        xnew <- object$x
    }
    if (is.null(Z) & (length(object$ind.drift) >0 )) {
        Z <- object$Tmatrix[, !object$ind.drift]
    }
    if (!is.null(ynew)) {
        coef.hold <- mKrig.coef(object, ynew,
                                collapseFixedEffect=collapseFixedEffect)
        c.coef <- coef.hold$c
        d.coef <- coef.hold$d
    }
    else {
        c.coef <- object$c
        d.coef <- object$d
    }
    # fixed part of the model this a polynomial of degree m-1
    # Tmatrix <- fields.mkpoly(xnew, m=object$m)
    # only do this if nt>0, i.e. there is a fixed part.
    #
    if( object$nt>0){
    if (derivative == 0) {
        if (drop.Z | object$nZ == 0) {
            # just evaluate polynomial and not the Z covariate
            temp1 <- fields.mkpoly(xnew, m = object$m) %*% 
              d.coef[object$ind.drift, ]
        }
        else {
          if( nrow( xnew) != nrow(as.matrix(Z)) ){
            stop(paste("number of rows of covariate Z",
                         nrow(as.matrix(Z)), 
             " is not the same as the number of locations",
                         nrow( xnew) )
                )
          }
            temp0 <-  cbind(fields.mkpoly(xnew, m = object$m),as.matrix(Z)) 
            temp1 <- temp0 %*% d.coef
        }
    }
    else {
        if (!drop.Z & object$nZ > 0) {
            stop("derivative not supported with Z covariate included")
        }
        temp1 <- fields.derivative.poly(xnew, m = object$m, d.coef[object$ind.drift, 
            ])
    }
    if (just.fixed) {
        return(temp1)
    }
    }  
    # add nonparametric part. Covariance basis functions
    # times coefficients.
    # syntax is the name of the function and then a list with
    # all the arguments. This allows for different covariance functions
    # that have been passed as their name.
    if (derivative == 0) {
        # argument list are the parameters and other options from mKrig
        #  locations and coefficients,
        temp2 <- do.call(object$cov.function.name, c(object$args, 
            list(x1 = xnew, x2 = object$knots, C = c.coef), cov.args))
    }
    else {
        temp2 <- do.call(object$cov.function.name, c(object$args, 
            list(x1 = xnew, x2 = object$knots, C = c.coef, derivative = derivative), 
            cov.args))
    }
    # add two parts together and coerce to vector
    if( object$nt>0){
      return((temp1 + temp2))
    }
    else{
      return(  temp2)
    }
}
