
# comparing covariance to the approx one if neagtive values set equal to zero


gridList<- list(  1:40,  1:40 )
gridPoints<- make.surface.grid( gridList)
# create weights in Fourier domain
obj<- circulantEmbeddingSetup( gridList, M=c(80,80),
                  cov.args=list( Covariance="Matern", theta=15, smoothness=1.0)
)



# this is just the exact covariance function that we started with. 
look<- fft( obj$wght, inverse=TRUE)
# note Imaginary part is machine zero
stats( c(Im( look)) )
# reduce to the grid and omit the imaginary part
look<- (array( Re(look), dim(obj$wght) ))[ 1:40, 1:40]
image.plot( look)
# as a function of distance 
x0<- cbind( 1,1)
D<- rdist( gridPoints, x0)
plot( D, c( look), pch=16, type="p", cex=.5)

#as an exact simulation will fail because some negative weights

stats( c( Re(obj$wght))) #  note some negative values

stats( c( Im(obj$wght))) # machine zero


wght0<-  obj$wght
# set the negative ones to zero
wght0[ Re(wght0)< 0 ]<- 0


#look at the covariance function for this approximation 

look0<- fft( wght0, inverse=TRUE)
look0<- (array( Re(look0), dim(obj$wght) ))[ 1:40, 1:40]
image.plot( look0) 
# not too bad!

# relative error
image.plot( (look- look0)/look )
# within 2%


# look at the covariance along one axis
matplot( cbind( look[,1], look0[,1]),
         type="l")
# some error at zero 

# adjust to have varaince of 1  we expect look0 to have a lower variance because of setting some 
# values to zero 

look0<- look0/ max( look0)

matplot( cbind( look[,1], look0[,1]),
         type="l")






