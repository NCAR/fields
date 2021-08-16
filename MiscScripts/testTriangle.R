setwd("~/Dropbox/Home/Repositories/fields")
source("image.plot.R")
source("colorBar.R")
source( "addColorBarTriangle.R")

x<- c(200, 300, 400)
y <- x + 100
z<-  outer( x, y, "+")
obj<- list( x=x, y=y, z=z)

dev.off()

image.plotNEW(obj , horizontal= TRUE, col=rainbow(9))

dev.off()
image.plotNEW(obj , horizontal= TRUE,col=topo.colors(9),
              lowerColor="grey10", upperColor="orange3"
              )
dev.off()
quartz( width=3,height=3)
image.plotNEW(obj , horizontal= FALSE,col=topo.colors(9),
              lowerColor="magenta", upperColor="orange3"
)

image.plotNEW(obj , horizontal= FALSE,
              LColor="magenta", RColor="magenta"
)