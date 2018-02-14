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
print.spatialProcess <- function(x, digits = 4, ...) {
  
  if (is.matrix(x$residuals)) {
    n <- nrow(x$residuals)
    NData <- ncol(x$residuals)
  }
  else {
    n <- length(x$residuals)
    NData <- 1
  }
  
  c1 <- "Number of Observations:"
  c2 <- n
  
  if (NData > 1) {
    c1 <- c(c1, "Number of data sets fit:")
    c2 <- c(c2, NData)
  }
  
  c1 <- c(c1, "Degree of polynomial in fixed part: ")
  
  
  if(x$m !=0 ){
    c2 <- c(c2, x$m - 1)
  }
  else{
    c2 <- c(c2, NA)
  }
  c1 <- c(c1, "Total number of parameters in fixed part: ")
  c2 <- c(c2, x$nt)
  if (x$nZ > 0) {
    c1 <- c(c1, "Number of additional covariates (Z)")
    c2 <- c(c2, x$nZ)
  }
  
  c1 <- c(c1, "MLE nugget variance ( sigma^2)")
  c2 <- c(c2, signif(x$sigma.MLE.FULL^2, digits))
    
  c1 <- c(c1, "MLE process variance (rho)")
  c2 <- c(c2, signif(x$rho.MLE.FULL, digits))
  
  
  c1 <- c(c1, "MLE range parameter (theta, units of distance): ")
  c2 <- c(c2, signif(x$theta.MLE, digits))
  
  c1 <- c(c1, "Approx 95% lower bound:  ")
  c2<-  c( c2, signif( x$theta.95CI[1], digits) )
  
  c1 <- c(c1, "           upper bound:  ")
  c2<- c( c2, signif( x$theta.95CI[2], digits) )
  
 
  if (!is.na(x$eff.df)) {
    c1 <- c(c1, " Approx.  degrees of freedom for curve")
    c2 <- c(c2, signif(x$eff.df, digits))
    if (length(x$trA.info) < x$np) {
      c1 <- c(c1, "   Standard Error of df estimate: ")
      c2 <- c(c2, signif(sd(x$trA.info)/sqrt(length(x$trA.info)), 
                         digits))
    }
  }
  c1 <- c(c1, "Nonzero entries in covariance")
  c2 <- c(c2, x$nonzero.entries)
  sum <- cbind(c1, c2)
  dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
  
  cat("Call:\n")
  dput(x$call)
  print(sum, quote = FALSE)
  cat(" ", fill = TRUE)
  cat(" Covariance Model:", x$cov.function, fill = TRUE)
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
  if (!is.na(x$eff.df)) {
    c1 <- c(c1, " Eff. degrees of freedom")
    c2 <- c(c2, signif(x$eff.df, digits))
    if (length(x$trA.info) < x$np) {
      c1 <- c(c1, "   Standard Error of estimate: ")
      c2 <- c(c2, signif(sd(x$trA.info)/sqrt(length(x$trA.info)), 
                         digits))
    }
  }
  
  invisible(x)
}
