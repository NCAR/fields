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
print.spatialProcessSummary <- function(x, digits = 4, ...) {
  
  cat("CALL:\n")
  dput(x$call)
  cat("\n")
  cat(" SUMMARY OF MODEL FIT:\n")
  print(x$summaryTable, quote = FALSE)
  cat("\n")
  cat(" ESTIMATED COEFFICIENTS FOR FIXED PART:", fill = TRUE)
  cat("\n")
  print(x$fixedEffectsTable, quote = FALSE)
  cat("\n")
  cat(" COVARIANCE MODEL:", x$cov.function, fill = TRUE)
  if (x$cov.function == "stationary.cov") {
    covName<- ifelse(is.null(x$args$Covariance), 
                     "Exponential", x$args$Covariance)
    cat("  Covariance function: ", covName, fill = TRUE)
  }
  if (!is.null(x$args)) {
    cat("   Non-default covariance arguments and their values ", 
        fill = TRUE)
    nlist <- as.character(names(x$args))
    NL <- length(nlist)
    for (k in 1:NL) {
      covItem<- x$args[[k]]
      cat( nlist[k], ":",
           fill = TRUE)
      if (object.size(covItem) > 10) {
        print(covItem)
      }
      else {
        cat("Too large to print value, size > 10 ...", 
            fill = TRUE)
      }
    }
  }
  cat("\n")
  cat(" SUMMARY FROM ML ESTIMATION:", fill=TRUE)
  if( !is.null( x$MLEpars)){
    cat("  parameters found from optim: ", x$MLEpars, fill=TRUE )
    cat(" tau and sigma found analytically", fill=TRUE)
    cat("\n")
    print( x$MLESummary)
  }
  else{
    cat(" lambda is fixed: ", fill=TRUE)
    cat(" tau and sigma  analytically  derived from lambda", fill=TRUE)
    cat("\n")
    print( x$MLESummary)
  }
 
  invisible(x)
}
