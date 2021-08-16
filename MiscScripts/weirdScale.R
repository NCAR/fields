



z<- outer( seq( 1,50,,30), seq( 1,50,,30) ,"+")



breaks<- round ( c(  1:9, 10^seq(1, 2,length.out=5) ), 1)
n<- length( breaks)
colTab<- tim.colors( n-1)

zInt<- matrix( cut( z, breaks, labels=FALSE),
              nrow(z), ncol(z)
              )

image.plot(zInt,
          col = colTab,
    axis.args = list(
             at = 1:n - .5,
         labels = format( round(breaks,1))
           )
)

    


