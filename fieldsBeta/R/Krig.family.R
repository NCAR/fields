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

Krig.check.xY <- function(x, Y, Z, weights, na.rm, 
    verbose = FALSE) {
    #
    # check for missing values in Y or X.
    #
    # save logical indicating where there are NA's
    # and check for NA's
    #
    ind <- is.na(Y)
    if (any(ind) & !na.rm) {
        stop("Need to remove missing values or use: na.rm=TRUE in the call")
    }
    #
    # coerce x to be a matrix
    x <- as.matrix(x)
    #
    # coerce Y to be a vector
    #
    Y <- as.matrix(Y)
    if (ncol(Y) != 1) {
        stop("Krig can not handle matrix Y data. See mKrig.")
    }
    #
    #default weights ( reciprocal variance of errors).
    #
    if (is.null(weights)) 
        weights <- rep(1, length(Y))
    #
    # check that dimensions agree
    #
    if (length(Y) != nrow(x)) {
        stop(" length of y and number of rows of x differ")
    }
    if (length(Y) != length(weights)) {
        stop(" length of y and weights differ")
    }
    #  if Z is not NULL coerce to be  a matrix
    # and check  # of rows
    if (verbose) {
        print(Z)
    }
    if (!is.null(Z)) {
        if (!is.matrix(Z)) {
            Z <- matrix(c(Z), ncol = 1)
        }
        if (length(Y) != nrow(Z)) {
            stop(" length of y and number of rows of Z differ")
        }
    }
    # if NAs can be removed then remove them and warn the user
    if (na.rm) {
        ind <- is.na(Y)
        if(all(ind)){
        	stop("Oops! All Y values are missing!")
        }
        if (any(ind)) {
            Y <- Y[!ind]
            x <- as.matrix(x[!ind, ])
            if (!is.null(Z)) {
                Z <- Z[!ind, ]
            }
            weights <- weights[!ind]
        }
    }
    #
    # check for NA's in x matrix -- there should not be any !
    if (any(c(is.na(x)))) {
        stop(" NA's in x matrix")
    }
    #
    # check for NA's in Z matrix
    if (!is.null(Z)) {
        if (any(c(is.na(Z)))) {
            stop(" NA's in Z matrix")
        }
    }
    #
    # verbose block
    if (verbose) {
        cat("Y:", fill = TRUE)
        print(Y)
        cat("x:", fill = TRUE)
        print(x)
        cat("weights:", fill = TRUE)
        cat(weights, fill = TRUE)
    }
    #
    # save x, weights  and Y w/o NAs
    N <- length(Y)
    return(list(N = N, y = Y, x = x, weights = weights, Z = Z, 
        NA.ind = ind))
}

"Krig.coef" <- function(out, lambda = out$lambda, 
    y = NULL, yM = NULL, verbose = FALSE) {
    #
    # NOTE default value of lambda used from Krig object.
    #
    # Determine whether to collapse onto means of replicates ( using y)
    # if the data has been passed use as the replicate means (yM) use that.
    # If both y and YM are null then just use out$yM
    # For readability of this function, all this tortured logic happens in
    #  Krig.ynew.
    #
    out2 <- Krig.ynew(out, y, yM)
    temp.yM <- out2$yM
    nt <- out$nt
    np <- out$np
    ndata <- ncol(temp.yM)
    u <- NA
    call.name <- out$cov.function.name
    if (verbose) {
        cat("dimension of yM in Krig.coef", fill = TRUE)
        print(dim(temp.yM))
    }
    #
    #   case when knots= unqiue x's
    # any lambda
    #
    if (out$decomp == "WBW") {
        # pad u with zeroes that corresond to null space basis functions
        # this makes it compatible with the DR decomposition.
        u <- rbind(matrix(0, nrow = out$nt, ncol = ndata), t(out$matrices$V) %*% 
            qr.q2ty(out$matrices$qr.T, out$W2 %d*% temp.yM))
        #
        #old code   beta <- out$matrices$G %*% ((1/(1 + lambda * out$matrices$D))%d*%u)
        #
        ind <- (nt + 1):np
        D2 <- out$matrices$D[ind]
        #
        # note use of efficient diagonal multiply in next line
        temp2 <- (D2/(1 + lambda * D2)) %d*% u[ind, ]
        beta2 <- out$matrices$V %*% temp2
        temp.c <- rbind(matrix(0, nrow = nt, ncol = ndata), beta2)
        temp.c <- qr.qy(out$matrices$qr.T, temp.c)
        temp.c <- out$W2 %d*% temp.c
        temp <- temp.yM - do.call(call.name, c(out$args, list(x1 = out$knots, 
            x2 = out$knots, C = temp.c)))
        temp <- out$W2 %d*% temp
        temp.d <- qr.coef(out$matrices$qr.T, temp)
    }
    
    if (out$decomp == "cholesky") {
        if (lambda != out$matrices$lambda) {
            stop("New lambda can not be used with cholesky decomposition")
        }
        Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
            list(x = out$knots, Z = out$ZM)))
        temp.d <- qr.coef(out$matrices$qr.VT, forwardsolve(out$matrices$Mc, 
            transpose = TRUE, temp.yM, upper.tri = TRUE))
        temp.c <- forwardsolve(out$matrices$Mc, transpose = TRUE, 
            temp.yM - Tmatrix %*% temp.d, upper.tri = TRUE)
        temp.c <- backsolve(out$matrices$Mc, temp.c)
    }
    
    return(list(c = temp.c, d = temp.d, tauHat.rep = out2$tauHat.rep, 
        tauHat.pure.error = out2$tauHat.pure.error, pure.ss = out2$pure.ss))
}

Krig.Amatrix <- function(object, x0 = object$x, lambda = NULL, 
    eval.correlation.model = FALSE, ...) {
    if (is.null(lambda)) {
        lambda <- object$lambda
    }
    M <- nrow(object$xM)
    N <- nrow(x0)
    # create output matrix
    out <- matrix(NA, N, M)
    #
    # loop through unique data locations predicting response
    # using unit vector
    # NOTE that the y vector has already been collapsed onto means.
    #
    for (k in 1:M) {
        ytemp <- rep(0, M)
        ytemp[k] <- 1
        out[, k] <- predict(object, x = x0, yM = ytemp, lambda = lambda, 
            eval.correlation.model = eval.correlation.model, 
            ...)
    }
    return(out)
}
"Krig.df.to.lambda" <- function(df, D, guess = 1, 
    tol = 1e-05) {
    if (is.list(D)) {
        D <- D$matrices$D
    }
    if (is.na(df)) 
        return(NA)
    if (df < sum(D == 0)) {
        warning("df too small to match with a lambda value")
        return(NA)
    }
    if (df > length(D)) {
        warning(" df too large to match a lambda value")
        return(NA)
    }
    l1 <- guess
    for (k in 1:25) {
        tr <- sum(1/(1 + l1 * D))
        if (tr <= df) 
            break
        l1 <- l1 * 4
    }
    l2 <- guess
    for (k in 1:25) {
        tr <- sum(1/(1 + l2 * D))
        if (tr >= df) 
            break
        l2 <- l2/4
    }
    info <- list(D = D, df = df, N = length(D))
    out <- bisection.search(log(l1), log(l2), Krig.fdf, tol = tol, 
        f.extra = info)$x
    +exp(out)
}

"Krig.engine.default" <- function(out, verbose = FALSE) {
    
    Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
        list(x = out$xM, Z = out$ZM)))
    if (verbose) {
        cat(" Model Matrix: spatial drift and Z", fill = TRUE)
        print(Tmatrix)
    }
    # Tmatrix premultiplied by sqrt of weights
    Tmatrix <- out$W2 %d*% Tmatrix
    qr.T <- qr(Tmatrix)
    if( qr.T$rank < ncol( Tmatrix)){
      stop("Regression matrix for fixed part of model is colinear")}
    #
    #verbose block
    if (verbose) {
        cat("first 5 rows of qr.T$qr", fill = TRUE)
        print(qr.T$qr[1:5, ])
    }
    #
    # find  Q_2 K Q_2^T  where K is the covariance matrix at the knot points
    #
    tempM <- t(out$W2 %d*% do.call(out$cov.function.name, c(out$args, 
        list(x1 = out$knots, x2 = out$knots))))
    tempM <- out$W2 %d*% tempM
    tempM <- qr.yq2(qr.T, tempM)
    tempM <- qr.q2ty(qr.T, tempM)
    np <- nrow(out$knots)
    nt <- (qr.T$rank)
    if (verbose) {
        cat("np, nt", np, nt, fill = TRUE)
    }
    #
    # Full set of decompositions for
    # estimator for nonzero lambda
    tempM <- eigen(tempM, symmetric = TRUE)
    D <- c(rep(0, nt), 1/tempM$values)
    #
    # verbose block
    if (verbose) {
        cat("eigen values:", fill = TRUE)
        print(D)
    }
    #
    # Find the transformed data vector used to
    # evaluate the solution, GCV, REML  at different lambdas
    #
    
    u <- c(rep(0, nt), t(tempM$vectors) %*% qr.q2ty(qr.T, c(out$W2 %d*% 
        out$yM)))
    if (verbose) {
        cat("u vector:", fill = TRUE)
        print(u)
    }
    #
    #
    return(list(D = D, qr.T = qr.T, decomp = "WBW", V = tempM$vectors, 
        u = u, nt = nt, np = np))
}

"Krig.engine.fixed" <- function(out, verbose = FALSE, 
    lambda = NA) {
    if (is.na(lambda)) 
        lambda <- out$lambda
    call.name <- out$cov.function.name
   
        Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
            list(x = out$knots, Z = out$ZM)))
        if (verbose) {
            cat("Tmatrix:", fill = TRUE)
            print(Tmatrix)
        }
        np <- nrow(out$knots)
        nt <- ncol(Tmatrix)
        # form K
        tempM <- do.call(call.name, c(out$args, list(x1 = out$knots, 
            x2 = out$knots)))
        # form M
        diag(tempM) <- (lambda/out$weightsM) + diag(tempM)
        #
        # find cholesky factor
        #  tempM = t(Mc)%*% Mc
        #  V=  Mc^{-T}
        # call cholesky but also add in the args supplied in Krig object.
        Mc <- do.call("chol", c(list(x = tempM), out$chol.args))
        VT <- forwardsolve(Mc, x = Tmatrix, transpose = TRUE, 
            upper.tri = TRUE)
        qr.VT <- qr(VT)
        # find GLS covariance matrix of null space parameters.
        Rinv <- solve(qr.R(qr.VT))
        Omega <- Rinv %*% t(Rinv)
        #
        # now do generalized least squares for d
        # and then find c.
        beta <- qr.coef(qr.VT, forwardsolve(Mc, transpose = TRUE, 
            out$yM, upper.tri = TRUE))
        if (verbose) {
          cat("beta fixed coefficients", fill=TRUE)
            print(beta)
        }
        c.coef <- forwardsolve(Mc, transpose = TRUE, out$yM - 
            Tmatrix %*% beta, upper.tri = TRUE)
        c.coef <- backsolve(Mc, c.coef)
        # return all the goodies,  include lambda as a check because
        # results are meaningless for other values of lambda
        return(list(qr.VT = qr.VT, d = c(beta), c = c(c.coef), 
            Mc = Mc, decomp = "cholesky", nt = nt, np = np, lambda.fixed = lambda, 
            Omega = Omega))
}

"Krig.fdf" <- function(llam, info) {
    sum(1/(1 + exp(llam) * info$D)) - info$df
}

"Krig.fgcv" <- function(lam, obj) {
    #
    # GCV that is leave-one-group out
    #
    lD <- obj$matrices$D * lam
    RSS <- sum(((obj$matrices$u * lD)/(1 + lD))^2)
    MSE <- RSS/length(lD)
    if ((obj$N - length(lD)) > 0) {
        MSE <- MSE + obj$pure.ss/(obj$N - length(lD))
    }
    trA <- sum(1/(1 + lD))
    den <- (1 - (obj$cost * (trA - obj$nt - obj$offset) + obj$nt)/length(lD))
    # If the denominator is negative then flag this as a bogus case
    # by making the GCV function 'infinity'
    #
    ifelse(den > 0, MSE/den^2, 1e20)
}

"Krig.fgcv.model" <- function(lam, obj) {
    lD <- obj$matrices$D * lam
    MSE <- sum(((obj$matrices$u * lD)/(1 + lD))^2)/length(lD)
    trA <- sum(1/(1 + lD))
    den <- (1 - (obj$cost * (trA - obj$nt - obj$offset) + obj$nt)/length(lD))
    ifelse(den > 0, obj$tauHat.pure.error^2 + MSE/den^2, 1e20)
}

"Krig.fgcv.one" <- function(lam, obj) {
    lD <- obj$matrices$D * lam
    RSS <- obj$pure.ss + sum(((obj$matrices$u * lD)/(1 + lD))^2)
    trA <- sum(1/(1 + lD))
    den <- 1 - (obj$cost * (trA - obj$nt - obj$offset) + obj$nt)/obj$N
    # If the denominator is negative then flag this as a bogus case
    # by making the GCV function 'infinity'
    #
    ifelse(den > 0, (RSS/obj$N)/den^2, 1e+20)
}

"Krig.flplike" <- function(lambda, obj) {
    #  - log profile likelihood for lambda
    # See section 3.4 from Nychka  Spatial Processes as Smoothers paper.
    # for equation and derivation
    D2 <- obj$matrices$D[obj$matrices$D > 0]
    u2 <- obj$matrices$u[obj$matrices$D > 0]
    lD <- D2 * lambda
    N2 <- length(D2)
    # MLE estimate of sigma for fixed lambda
    sigma.MLE <- (sum((D2 * (u2)^2)/(1 + lD)))/N2
    #
    # ln determinant of    K + lambda*WI
    lnDetCov <- -sum(log(D2/(1 + lD)))
    
    -1 * (-N2/2 - log(2 * pi) * (N2/2) - (N2/2) * log(sigma.MLE) - 
        (1/2) * lnDetCov)
      
    
}

"Krig.fs2hat" <- function(lam, obj) {
    lD  <- obj$matrices$D * lam
    RSS <- obj$pure.ss + sum(((obj$matrices$u * lD)/(1 + lD))^2)
    den <- obj$N - (sum(1/(1 + lD)) + obj$offset)
    if (den < 0) {
        return(NA)
    }
    else {
        RSS/(den)
    }
}

"Krig.ftrace" <- function(lam, D) {
    sum(1/(1 + lam * D))
}

"Krig.make.W" <- function(out, verbose = FALSE) {
    if (verbose) {
        cat("W", fill = TRUE)
        print(out$W)
    }
    if (out$nondiag.W) {
        #
        # create W from scratch or grab it from passed object
        if (is.null(out$W)) {
            if (verbose) {
                print(out$wght.function.name)
            }
            W <- do.call(out$wght.function.name, c(list(x = out$xM), 
                out$wght.args))
            #       adjust W based on diagonal weight terms
            #
            W <- sqrt(out$weightsM) * t(sqrt(out$weightsM) * 
                W)
        }
        else {
            W <- out$W
        }
        #
        # symmetric square root
        temp <- eigen(W, symmetric = TRUE)
        W2 <- temp$vectors %*% diag(sqrt(temp$values)) %*% t(temp$vectors)
        return(list(W = W, W2 = W2))
    }
    else {
        #
        #  These are created only for use with default method to stay
        #   consistent with nondiagonal elements.
        if (out$fixed.model) {
            return(list(W = NULL, W2 = NULL))
        }
        else {
            return(list(W = out$weightsM, W2 = sqrt(out$weightsM)))
        }
    }
}

"Krig.make.Wi" <- function(out, verbose = FALSE) {
    #
    # If a weight matrix has been passed use it.
    #
    # Note that in either case the weight matrix assumes that
    # replicate observations have been collapses to the means.
    #
    if (out$nondiag.W) {
        temp <- eigen(out$W, symmetric = TRUE)
        Wi <- temp$vectors %*% diag(1/(temp$values)) %*% t(temp$vectors)
        W2i <- temp$vectors %*% diag(1/sqrt(temp$values)) %*% 
            t(temp$vectors)
        return(list(Wi = Wi, W2i = W2i))
    }
    else {
        #
        #  These are created only for use with default method to stay
        # consistent with nondiagonal elements.
        return(list(Wi = 1/out$weightsM, W2i = 1/sqrt(out$weightsM)))
    }
}

"Krig.make.u" <- function(out, y = NULL, yM = NULL, 
    verbose = FALSE) {
    #
    # Determine whether to collapse onto means of replicates ( using y)
    # if the data has been passed use as the replicate means (yM) use that.
    # If both y and YM are null then just use out$yM
    # For readability of this function, all this tortured logic happens in
    #  Krig.ynew.
    #
    out2 <- Krig.ynew(out, y, yM)
    temp.yM <- out2$yM
    nt <- out$nt
    np <- out$np
    ndata <- ncol(temp.yM)
    u <- NA
    call.name <- out$cov.function.name
    if (verbose) {
        cat("dimension of yM in Krig.coef", fill = TRUE)
        print(dim(temp.yM))
    }
    #
    #   case when knots= unqiue x's
    # any lambda
    #
        u <- rbind(matrix(0, nrow = out$nt, ncol = ndata), t(out$matrices$V) %*% 
            qr.q2ty(out$matrices$qr.T, out$W2 %d*% temp.yM))
    
    return(list(u = u, tauHat.rep = out2$tauHat.rep, tauHat.pure.error = out2$tauHat.pure.error, 
        pure.ss = out2$pure.ss))
}

Krig.null.function <- function(x, Z = NULL, drop.Z = FALSE, 
    m) {
    # default function to create matrix for fixed part of model
    #  x, Z, and drop.Z are required
    #  Note that the degree of the polynomial is by convention (m-1)
    # returned matrix must have the columns from Z last!
    #
    if (is.null(Z) | drop.Z) {
        return(fields.mkpoly(x, m = m))
    }
    else {
        return(cbind(fields.mkpoly(x, m = m), Z))
    }
}

Krig.parameters <- function(obj, mle.calc = obj$mle.calc) {
    # if nondiag W is supplied then use it.
    # otherwise assume a diagonal set of weights.
    #
    # NOTE: calculation of  tauHat involves full set of obs
    # not those colllapsed to the mean.
    if (obj$nondiag.W) {
        tauHat.GCV <- sqrt(sum((obj$W2 %d*% obj$residuals)^2)/(length(obj$y) - 
            obj$eff.df))
    }
    else {
        tauHat.GCV <- sqrt(sum((obj$weights * obj$residuals^2)/(length(obj$y) - 
            obj$eff.df)))
    }
    if (mle.calc) {
        sigma.MLE <- sum(c(obj$c) * c(obj$yM))/obj$N
        # set sigma estimate to zero if negtive. Typically this
        # is an issue of machine precision and very small negative value.
        sigma.MLE <- ifelse(sigma.MLE < 0, 0, sigma.MLE)
        
        #    commented out code for debugging ...
        #      if( sigma.MLE< 0) {
        #        stop('problems computing sigma.MLE')}
        # commented out is the REML estimate -- lose null space df because of
        # the restiction to orthogonal subspace of T.
        # sigmahat<- sigma.MLE <- sum(obj$c * obj$yM)/(obj$N - obj$nt)
        # .
        sigmahat <- sigma.MLE
        tauHat.MLE <- sqrt(sigma.MLE * obj$lambda)
    }
    else {
        sigmahat <- sigma.MLE <- tauHat.MLE <- NA
    }
    list(tauHat.GCV = tauHat.GCV, sigma.MLE = sigma.MLE, tauHat.MLE = tauHat.MLE, 
        sigmahat = sigmahat)
}

"Krig.replicates" <- function(out=NULL, x,y, Z=NULL, weights=rep( 1, length(y)),
                               digits=8,
                               verbose = FALSE) {
    if( is.null(out)){
      out<- list( x=x, y=y, N= length(y), Z=Z, weights=weights)
    }
    rep.info <- cat.matrix(out$x, digits=digits)
    if (verbose) {
        cat("replication info", fill = TRUE)
        print(rep.info)
    }
    # If no replicates are found then reset output list to reflect this condition
    uniquerows <- !duplicated(rep.info)
    if (sum(uniquerows) == out$N) {
        tauHat.rep <- NA
        tauHat.pure.error <- NA
        pure.ss <- 0
        # coerce 'y' data vector as a single column matrix
        yM <- as.matrix(out$y)
        weightsM <- out$weights
        xM <- as.matrix(out$x[uniquerows, ])
        # coerce ZM to matrix
        if (!is.null(out$Z)) {
            ZM <- as.matrix(out$Z)
        }
        else {
            ZM <- NULL
        }
    }
    # collapse over spatial replicates
    else {
        rep.info.aov <- fast.1way(rep.info, out$y, out$weights)
        tauHat.pure.error <- sqrt(rep.info.aov$MSE)
        tauHat.rep <- tauHat.pure.error
        # copy  replicate means as a single column matrix
        yM <- as.matrix(rep.info.aov$means)
        weightsM <- rep.info.aov$w.means
        xM <- as.matrix(out$x[uniquerows, ])
        # choose some Z's for replicate group means
        if (!is.null(out$Z)) {
            ZM <- as.matrix(out$Z[uniquerows, ])
        }
        else {
            ZM <- NULL
        }
        pure.ss <- rep.info.aov$SSE
        if (verbose) 
            print(rep.info.aov)
    }
    return(list(yM = yM, xM = xM, ZM = ZM, weightsM = weightsM, 
        uniquerows = uniquerows, tauHat.rep = tauHat.rep, tauHat.pure.error = tauHat.pure.error, 
        pure.ss = pure.ss, rep.info = rep.info))
}

Krig.transform.xY <- function(obj, knots=NA,  verbose = FALSE) {
    # find all replcates and  collapse to unique locations and mean response
    # and pooled variances and weights.
    out <- Krig.replicates(obj, verbose = verbose)
    if (verbose) {
        cat("yM from Krig.transform.xY", fill = TRUE)
        print(out$yM)
    }
    #
    # save information about knots.
        out$knots <- out$xM
        out$mle.calc <- TRUE
        out$knot.model <- FALSE
    #
    # scale x, knot locations and  save transformation info
    #
    out$xM <- transformx(out$xM, obj$scale.type, obj$x.center, 
        obj$x.scale)
    out$transform <- attributes(out$xM)
    out$knots <- scale(out$knots, center = out$transform$x.center, 
        scale = out$transform$x.scale)
    #
    #
    #verbose block
    #
    if (verbose) {
        cat("transform", fill = TRUE)
        print(out$transform)
    }
    if (verbose) {
        cat("knots in transformed scale", fill = TRUE)
        print(knots)
    }
    return(out)
}

"Krig.updateY" <- function(out, Y, verbose = FALSE, 
    yM = NA) {
    #
    #
    if (is.na(yM[1])) {
        out2 <- Krig.ynew(out, Y)
    }
    else {
        out2 <- list(yM = yM, tauHat.rep = NA, tauHat.pure.error = NA, 
            pure.ss = NA)
    }
    if (verbose) {
        print(out2)
    }
    #
    # Note how matrices are grabbed from the Krig object
    #
    if (verbose){ 
        cat("Type of decomposition", out$decomp, fill = TRUE)
     }
        #### decomposition of Q2TKQ2
        u <- c(rep(0, out$nt), t(out$matrices$V) %*% qr.q2ty(out$matrices$qr.T, 
            out$W2 %d*% out2$yM))
        if (verbose) 
            cat("u", u, fill = TRUE)
        #
        # pure error in this case from 1way ANOVA
        #
        if (verbose) {
            cat("pure.ss", fill = TRUE)
            print(out2$pure.ss)
        }

    out2$u <- u
    out2
}
Krig.which.lambda <- function(out) {
    #
    # determine the method for finding lambda
    #  Note order
    # default is to do 'gcv/REML'
    out2 <- list()
    # copy all all parameters to out2 just to make this
    # easier to read.
    out2$method <- out$method
    out2$lambda.est <- NA
    out2$lambda <- out$lambda
    out2$eff.df <- out$eff.df
    out2$sigma <- out$sigma
    out2$tau2 <- out$tau2
    if (!is.na(out2$lambda) | !is.na(out2$eff.df)) {
        #
        # this indicates lambda has been supplied and leads to
        # the cholesky type computational approaches
        #        -- but only if GCV is FALSE
        #
        out2$method <- "user"
    }
    out2$GCV <- out$GCV
    if (!is.na(out2$eff.df)) {
        #
        # this indicates df has been supplied and needs
        # GCV to be true to compute the lambda
        # that matches the df
        #
        out2$GCV <- TRUE
    }
    if (!is.na(out2$sigma) & !is.na(out2$tau2)) {
        out2$method <- "user"
        out2$lambda <- out2$tau2/out2$sigma
    }
    #
    # NOTE: method='user' means that a value of lambda has been supplied
    #        and so GCV etc to determine lambda is not needed.
    #  gcv TRUE means that the decompositions will be done to
    #    evaluate the estimate at arbitrary lambda (and also be
    #    able to compute the effective degrees of freedom).
    #
    #    The fixed lambda calculations are very efficient but
    #    do not make it feasible for GCV/REML  or effective degrees of
    #    freedom calculations.
    #
    out2$fixed.model <- (out2$method == "user") & (!out2$GCV)
    #
    return(out2)
}

"Krig.ynew" <- function(out, y = NULL, yM = NULL) {
    #
    # calculates the collapsed y (weighted) mean vector based on the
    # X matrix and weights from the out object.
    # or just passes through the collapsed mean data if passed.
    #
    #
    # If there are no replicated obs. then return the full vector
    # pure error ss is zero
    #
    tauHat.rep <- NA
    tauHat.pure.error <- NA
    pure.ss <- 0
    # if no y's are given then it is assumed that one should use the
    # yM from the original data used to create the Krig object
    if (is.null(yM) & is.null(y)) {
        yM <- out$yM
    }
    #
    # case when yM is passed no calculations are needed
    #
    if (!is.null(yM)) {
        return(list(yM = as.matrix(yM), tauHat.rep = NA, tauHat.pure.error = NA, 
            pure.ss = 0))
    }
    #
    # no reps case
    #
    if (length(unique(out$rep.info)) == out$N) {
        return(list(yM = as.matrix(y), tauHat.rep = NA, tauHat.pure.error = NA, 
            pure.ss = 0))
    }
    #
    #  check that y is the right length
    #
    if (length(y) != out$N) {
        stop(" the new y vector is the wrong length!")
    }
    #
    # case when full y data is passed and replicate means need to be found
    #
    if (length(unique(out$rep.info)) < out$N) {
        #
        # calculate means by pooling Replicated observations but use the
        # the right weighting.
        #
        rep.info.aov <- fast.1way(out$rep.info, y, out$weights)[c("means", 
            "MSE", "SSE")]
        tauHat.pure.error <- sqrt(rep.info.aov$MSE)
        tauHat.rep <- tauHat.pure.error
        return(list(yM = rep.info.aov$means, tauHat.rep = tauHat.rep, 
            tauHat.pure.error = tauHat.pure.error, pure.ss = rep.info.aov$SSE))
    }
}
