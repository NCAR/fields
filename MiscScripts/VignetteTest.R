library( fields)
source("fields12.3Patch.R")

data( COmonthlyMet)
# predicting average daily minimum temps for spring in Colorado
obj<- spatialProcess( CO.loc, CO.tmin.MAM.climate, Z= CO.elev)

out.p<-predictSurface( obj, grid.list=CO.Grid,
                             ZGrid= CO.elevGrid, extrap=TRUE)

image.plot( out.p, col=larry.colors())
US(add=TRUE, col="grey")
contour( CO.elevGrid, add=TRUE, levels=seq(1000,3000,,5), col="black")
title("Average Spring daily min. temp in CO")

out.p<-predictSurfaceSE( obj, grid.list=CO.Grid,
                         ZGrid= CO.elevGrid,
                         extrap=TRUE, verbose=FALSE) 

image.plot( out.p, col=larry.colors())
US(add=TRUE, col="grey")

data( COmonthlyMet)

# predicting average daily minimum temps for spring in Colorado
obj<- spatialProcess( CO.loc, CO.tmin.MAM.climate, Z= CO.elev) 
out.p<-predictSurface( obj, grid.list=CO.Grid, 
                            ZGrid= CO.elevGrid, extrap=TRUE)
image.plot( out.p, col=larry.colors())

US(add=TRUE, col="grey")
contour( CO.elevGrid, add=TRUE, levels=seq(1000,3000,,5), col="black") 
title("Average Spring daily min. temp in CO")

out.p<-predictSurfaceSE.default( obj, 
                                 grid.list=CO.Grid,
                                 ZGrid= CO.elevGrid, 
                                 extrap=TRUE)

#ZGrid= # error drop.Z not supported ??
predictSE.mKrig( obj)
predictSE.mKrig( obj, obj$x[1:10,], Z= CO.elev[1:10], verbose=TRUE)
                 
out.p2<-predictSE.mKrig( obj, make.surface.grid(CO.Grid), 
                  Z= as.matrix(c(CO.elevGrid$z) )
                  

image.plot( out.p, col=larry.colors())

fit<- spatialProcess(ChicagoO3$x,ChicagoO3$y,
                     aRange=1, lambda=.1)
predictSE.mKrig( fit,ChicagoO3$x )
