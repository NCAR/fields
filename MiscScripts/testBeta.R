library( fields)


data( ozone2)
# x is a two column matrix where each row is a location in lon/lat 
# coordinates
s<- ozone2$lon.lat
# y is a vector of ozone measurements at day 16. Note some missing values. 
z<- ozone2$y[16,]
good<- !is.na( z)
s<- s[good,]
z<- z[good]

############################################
#######  developing spatialProcess


obj0 <- spatialProcess(s,z,
                        aRange=.8, lambda=.1)
obj00 <- spatialProcess(s,z,
        cov.args=list(aRange=.8, lambda=.1)
)

obj2<- spatialProcess(s,z,
     cov.params.start = list( lambda= .1),
                      aRange=.8 )

obj3<- spatialProcess(s,z,
      cov.params.start = list( lambda= .1, aRange=.8))
 
rbind( obj0$summary,obj2$summary, obj3$summary )                       

parGrid<- data.frame( lambda= 10**seq( -3, 1, length.out=10) )
obj4<- spatialProcess(s,z, parGrid=parGrid, aRange=.8)
obj4A<- spatialProcess(s,z,  aRange=.8)

rbind( obj4$summary, obj4A$summary)

parGrid<- expand.grid( lambda= 10**seq( -3, 1, length.out=15) ,
                       aRange=seq( .5,1.0, length.out=15)
                       )
obj5<- spatialProcess(s,z, parGrid=parGrid)

obj6<- spatialProcess(s,z)

lambdaHat<- obj6$MLEInfo$pars.MLE[1]
aHat<- obj6$MLEInfo$pars.MLE[2]

obj7<- spatialProcess(s,z, aRange=aHat,
          cov.params.start= list(lambda = lambdaHat),
                                             verbose=FALSE)

look<- profileMLE( obj6, "aRange", 
        cov.params.start = list( lambda=lambdaHat),
          parGrid= data.frame( aRange=aHat),
      verbose=FALSE)

look2<- profileMLE( obj6, "aRange", 
                   cov.params.start = list( lambda=lambdaHat),
                   parGrid= data.frame( aRange=aHat),
                   verbose=FALSE, fast=TRUE)

obj7$summary
look$summary

obj8<- spatialProcess(s,z, profileLambda = TRUE, profileARange = TRUE)
obj9<- spatialProcess(s,z, profileLambda = TRUE, profileARange = TRUE)
obj10 <-profileMLE(obj9,"aRange", fast=TRUE)





