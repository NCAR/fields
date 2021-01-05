library( fields)


data( ozone2)
# x is a two column matrix where each row is a location in lon/lat 
# coordinates
s<- ozone2$lon.lat
# y is a vector of ozone measurements at day 16. Note some missing values. 
z<- ozone2$y[16,]
good<- !is.na( y)
s<- s[good,]
z<- z[good]

############################################
#######  developing spatialProcess


obj0 <- spatialProcess(s,z,
                        aRange=.8, lambda=.1)
look<- summary.spatialProcess(obj0)
plot.spatialProcess( obj0)

obj00 <- spatialProcess(s,z,
                       cov.args= list( aRange=.8, lambda=.1)
                       ) 

obj2<- spatialProcess(s,z,  aRange=.8 )
plot.spatialProcess( obj2)

obj2B<- spatialProcess(s,z,  aRange=.8,profileLambda = TRUE )
plot.spatialProcess( obj2B)

obj2C<- spatialProcess(s,z,  aRange=.8, profileLambda = TRUE, 
                       gridLambda= 10**seq( -3,0,,50) ,mKrig.args= list( iseed=34,
                                                                         NtrA=200)
                       )
plot.spatialProcess( obj2C)


obj3<- spatialProcess(s,z)
plot.spatialProcess( obj3)

obj4 <- spatialProcess(s,z, profileLambda=TRUE)

objGCV <- spatialProcess(s,z, aRange=1.0, gridLambda = seq(.1,1,,10 ),
                         GCV=TRUE, verbose=TRUE)

objGCV <- spatialProcess(s,z, aRange=1.0, gridLambda = seq(.1,1,,10 ),
                         verbose=TRUE, profileLambda = TRUE)


obj3B<- spatialProcess(s,z, 
                       cov.params.start= list( aRange=.8, lambda=.1),
                       profileARange = TRUE, gridN= 20
                       )

plot.spatialProcess( obj3B)

obj4<- spatialProcess(s, z, profileLambda = TRUE)
plot.spatialProcess( obj4 )

fit<- spatialProcess(s,z)

summary( fit)

predict( fit)
predictSE( fit)

look1<- predictSurface(fit)
image.plot( look1)
look2<- predictSurfaceSE( fit)
image.plot( look2)


data( ozone2)
# x is a two column matrix where each row is a location in lon/lat 
# coordinates
s<- ozone2$lon.lat
# y is a vector of ozone measurements at day 16. Note some missing values. 
z<- ozone2$y[16,]
good<- !is.na( z)
s<- s[good,]
z<- z[good]


objTest<- fastTps( s,z, aRange= 1.0, profileLambda=TRUE)
set.panel(2,2)
plot(objTest)
set.panel()
surface( objTest)
look<- sim.fastTps.approx(objTest, M=5)
look<- sim.spatialProcess(objTest, s, M=5)

look<- sim.mKrig.approx(objTest, s, M=5)



