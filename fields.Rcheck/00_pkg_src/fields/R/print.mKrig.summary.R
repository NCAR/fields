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