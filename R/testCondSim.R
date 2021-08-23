
setwd("~/Dropbox/Home/Repositories/fields")
library( fields)
data( "ozone2")
mKrigObject<- spatialProcess( ozone2$lon.lat, ozone2$y[16,],
                              smoothness=.5)

nx<- 80
ny<- 80


xGridList<- fields.x.to.grid( mKrigObject$x, nx=nx, ny=ny)
xGrid<- make.surface.grid( xGridList)
allTime0<- system.time(
 look0<- sim.spatialProcess(mKrigObject, xp= xGrid, M=10)
)

print( allTime0)


allTime<- system.time(
  look<- simLocal.mKrig(mKrigObject, M=10,nx=nx, ny=ny,
                     gridRefinement = 2,
                     np=3)
) 

print( allTime)
print( look$timing)


allTime<- system.time(
  look<- simLocal.mKrig(mKrigObject, M=10,nx=64*3, ny=64*3,
                     gridRefinement = 1,
                     np=3)
) 
print( allTime)
print( look$timing)
user  system elapsed 
# 31.550   1.677  33.058 
# > print( look$timing)
# user.self sys.self elapsed user.child sys.child
# timeSetup1     2.040    0.134   2.174          0         0
# timeSetup2     0.002    0.002   0.004          0         0
# timePred1      0.171    0.040   0.200          0         0
# timePred2      0.168    0.027   0.181          0         0
# 



n<- 2500
s<- matrix( runif(n*2 ),n,2)
y<-  rnorm( n)

mKrigObject1<- spatialProcess( s,y,aRange=.1, lambda=.1,
                               smoothness=.5
)

allTime1<- system.time(
  look1<- simLocal.mKrig(mKrigObject1, M=10,nx=128, ny=128,
                     gridRefinement = 3,
                     np=3)
)

print( allTime1)
print( look1$timing)


mKrigObject2<- fastTps( s,y, lambda=.1, theta=.25)

allTime2<- system.time(
  look2<- simLocal.mKrig(mKrigObject2, M=10,nx=128*3, ny=128*3,
                      gridRefinement = 1, delta=.25,
                      np=3)
)

print( allTime2)
print( look2$timing)
# > print( allTime2)
# user  system elapsed 
# 26.132   5.833  32.833 
# > print( look2$timing)
# user.self sys.self elapsed user.child sys.child
# timeSetup1     5.836    4.964  11.635          0         0
# timeSetup2     0.031    0.013   0.048          0         0
# timePred1      1.275    0.013   1.288          0         0
# timePred2      1.316    0.009   1.324          0         0
# > 
