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
"simNN.mKrig" <- function(mKrigObject,  
    predictionGridList = NULL,
    simulationGridList = NULL, 
        gridRefinement = 1, 
                    np = 2, 
                     M = 1,
                    nx = 256,
                    ny = 256, 
               verbose = FALSE,
                          ...)
    {
   if (ncol(mKrigObject$x) != 2) {
        stop("conditional simulation only implemented for 2 dimensions")
    }
    nObs<- nrow( mKrigObject$x)
    # create prediction set of points based on what is passed
    if (is.null(predictionGridList)) {
        # these adjustments insure there are enough grid
        # points beyond the range of the locations. 
        xr<- range(mKrigObject$x[,1] )
        dx<- (xr[2]-xr[1])/(nx- 1 - 2*np )
        yr<- range(mKrigObject$x[,2] )
        dy<- (yr[2]-yr[1])/(ny - 1 - 2*np)
        predictionGridList<- list( x = seq( xr[1] - dx*np, xr[2] +dx*np, 
                                          length.out=nx), 
                                   y = seq( yr[1]- dy*np, yr[2] +dy*np,
                                          length.out=ny)
                                 )
    }
    else{
        nx<- length( predictionGridList$x)
        ny<- length( predictionGridList$y)
    }
    
    #
     if (is.null(simulationGridList)) {
         simulationGridList<- list( x= seq( min(predictionGridList$x), 
                                           max(predictionGridList$x),
                                           length.out=nx*gridRefinement ),
                                    y= seq( min(predictionGridList$y), 
                                            max(predictionGridList$y),
                                            length.out=ny*gridRefinement )
         )
         indexSubset<-  list( x= match(predictionGridList$x,
                                       simulationGridList$x),
                              y = match(predictionGridList$y,
                                        simulationGridList$y))
         if( any( is.na( indexSubset))){
             stop("predict grid is not a subset 
                  of the simulation grid")
         }
     }
    
    tau <-    mKrigObject$summary["tau"]
    sigma2 <- mKrigObject$summary["sigma2"]
    aRange<-  mKrigObject$summary["aRange"]
    Covariance <- mKrigObject$args$Covariance
    # wipe out some extraneous components that are not used by the Covariance
    # function.
    covArgs0 <- mKrigObject$args
    covArgs0$Covariance<- NULL
    covArgs0$distMat <- NULL
    covArgs0$onlyUpper<- NULL
    covArgs0$aRange<- NULL
    
    
    #
    # set up various sizes of arrays
    nObs <- nrow(mKrigObject$x)
    if (verbose) {       
        cat("nObs, tau, sigma2", nObs, tau, sigma2, fill = TRUE)
    }
    timeSetup1<- system.time(
    # set up object for simulating on a grid using circulant embedding
    CEObject<- circulantEmbeddingSetup(simulationGridList,
                                   cov.function = mKrigObject$cov.function,
                                       cov.args = mKrigObject$args )
    )
    if (verbose) {
        cat("dim of full circulant matrix ", dim(CEObject$wght), 
            fill = TRUE)
    }
    timeSetup2<- system.time(
    offGridObject<- offGridWeights( mKrigObject$x,
                                    simulationGridList,
                                    mKrigObject,
                                    np=np
                                    )
    )
    #
    # find conditional mean field from initial fit
    timePred1<- system.time(
      hHat <- predictSurface(mKrigObject,
                           grid.list = predictionGridList,
                            ...)$z
    )
    # setup output array to hold ensemble
    out <- array(NA, c( nx, ny, M))
    # empty image object to hold simulated fields
    
    ##########################################################################################
    ### begin the big loop
    ##########################################################################################
    sdNugget<- tau* sqrt(1/mKrigObject$weights)
    sdPredictionError<- sqrt(offGridObject$predictionVariance)
     for (k in 1:M) {
        cat(k, " ")
        # simulate full field
        
        hTrue<- sqrt(sigma2) * circulantEmbedding(CEObject)
        #
        # NOTE: fixed part of model (null space) does not need to be simulated
        # because the  estimator is unbiased for this part.
        # the variability is still captured because the fixed part
        # is still estimated as part of the predict step below
        #
        #   bilinear interpolation to approximate values at data locations
        #
        hData <- offGridObject$B%*%c(hTrue) + 
                   sdPredictionError*rnorm(nObs) 
        ySynthetic <- hData + sdNugget*rnorm(nObs)
        if (verbose) {
            cat("stats for synthetic values", fill = TRUE)
            print(t(stats(ySynthetic)))
        }
        # predict at grid using these data
        # and subtract from synthetic 'true' value
        #
        timePred2<-system.time(
        spatialError <- predictSurface(
                               mKrigObject,  
                   grid.list = predictionGridList, 
                        ynew = ySynthetic, 
                               ...)$z
        )
        spatialError <- spatialError - hTrue[ indexSubset$x,indexSubset$y] 
        # add the error to the actual estimate  (conditional mean)
        out[,, k] <- hHat + spatialError
    }
    return(list(x = predictionGridList$x,
                y = predictionGridList$y,
                z = out, 
                timing= rbind( timeSetup1, timeSetup2,
                               timePred1,timePred2),
             call = match.call()))
}

