library( fields)

#  colorado climate data 
  data(COmonthlyMet)
 
  x<- CO.loc
  y<- CO.tmin.MAM.climate
  elev<- CO.elev
  good<- !is.na( y)
  x<- x[good,]
  y<- y[good]
  elev<- elev[good]
 
    NGRID <- 50 
 # get elevations on a grid (will use these later) 
 #   lonRange<- range(x[,1])
 #   latRange<- range(x[,2])
 #  COGrid<- list(x=seq(lonRange[1], lonRange[2],length.out=NGRID ),
 #                y=seq(latRange[1], latRange[2],length.out=NGRID )
 #  )
    
   COGrid<- list(x = seq( -109, -101,,NGRID), 
                  y = seq(36.5, 41.5,,NGRID))
   COGridPoints<- make.surface.grid( COGrid )
   data( RMelevation)
   COElevGrid<- interp.surface.grid( RMelevation, COGrid )
#  Native Res:   COElevGrid<- crop.image( RMelevation, x )
 
  
# take a look at the data
  quilt.plot( x, y) 
  US( add=TRUE)
  
  plot( elev, y)
  plot( x[,2], y)
  
  X<- cbind( x, elev) 
# 
 lmObj<- lm( y ~ lon+lat +elev, data=X )  
 summary( lmObj)
  quilt.plot( x, lmObj$residuals) 
  US( add=TRUE)

   fit0E<- Tps( x,y, Z=elev- mean(elev))
   
   
# fit a Kriging estimator Matern covariance smoothness =1.0
# range and nugget estimated by maxmimum likelihood
  fit1<- spatialProcess( x,y)
  set.panel(2,2)
  plot( fit1)
  
  surface( fit1)
  US( add=TRUE)
  
  # take a look at residuals
  plot( elev, fit1$residuals)
  lmObj<- lm( fit1$residuals ~ elev)
  abline( lmObj, col="red", lwd=3)
  
  
  
# with elevations  
   fit1E<- spatialProcess( x,y, Z= elev)
# the clean way for fields >=9.0   
   sur0<- predictSurface( fit1E, grid.list=COGrid, ZGrid=COElevGrid )
# use this for fields <= 9.0 
   sur0<- predict( fit1E, x=COGridPoints, Z=c(COElevGrid$z) ) 
   sur0<- as.surface( COGrid, sur0)
#   
   image.plot( sur0, col=terrain.colors(256))
   sur0Smooth<- predictSurface( fit1E, grid.list=COGrid,  drop.Z= TRUE)
   image.plot( sur0Smooth)
   
# uncertainty 40 draws from posterior distribution
  set.seed(123)
  # next commmand takes a minute or so
  SEout<- sim.spatialProcess( fit1E, xp=COGridPoints,
                              Z= c(COElevGrid$z), M= 40, drop.Z=TRUE)

 
  set.panel( 3,3)                            
  image.plot( sur0Smooth, 
              zlim=c(3.5,14.5), col=tim.colors(256))
  contour( sur0Smooth, levels= 9, lwd=3, col="grey", add=TRUE) 
  par( mar=c(3,3,1,1))  
  for( k in 1:8){
      image( as.surface( COGrid, SEout[k,]), 
               zlim =c(3.5,14.5),col=tim.colors(256), axes=FALSE)
      contour( as.surface( COGrid, SEout[k,]),
               lwd=3, col="grey",level=9, add=TRUE)
               title( k, adj=0, cex=2)
  }
  
    
   surSE<- apply( SEout, 2, sd )
 
   image.plot( as.surface( COGrid, surSE))
   points( x, col="magenta", pch=16) 
   contour( COElevGrid, 
                level= c(2000, 3000), add= TRUE, col="grey30", lwd=3)
   US( add=TRUE, col="grey", lwd=2)
   

   
  
  
  
 
