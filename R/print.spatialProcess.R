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
  
  c1 <- c(c1, "Degree of polynomial null space ( base model):")
  
  
  if(x$m !=0 ){
    c2 <- c(c2, x$m - 1)
  }
  else{
    c2 <- c(c2, NA)
  }
  c1 <- c(c1, "Total number of parameters in base model")
  c2 <- c(c2, x$nt)
  if (x$nZ > 0) {
    c1 <- c(c1, "Number of additional covariates (Z)")
    c2 <- c(c2, x$nZ)
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
  c1 <- c(c1, "Smoothing parameter")
  c2 <- c(c2, signif(x$lambda.fixed, digits))
  
  c1 <- c(c1, "MLE sigma")
  c2 <- c(c2, signif(x$sigma.MLE.FULL, digits))
    
  c1 <- c(c1, "MLE rho")
  c2 <- c(c2, signif(x$rho.MLE.FULL, digits))
  
  c1 <- c(c1, "MLE lambda = MLE sigma^2 / MLE rho")
  c2 <- c(c2, signif(x$lambda.MLE, digits))
  
  c1 <- c(c1, "MLE theta")
  c2 <- c(c2, signif(x$theta.MLE, digits))
  
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
  invisible(x)
}
