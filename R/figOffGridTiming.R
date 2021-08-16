setwd("~/Dropbox/Home/Repositories/fields")
rm(list=ls())

library(fields)
library( tictoc)
source("timeOffGrid.R")
source("timeCEWithN.R")


np<- 4
out<- NULL
for(  m in c( 128,  256, 512,  2048)){
  for( n in c( 400,  1600, 6400)){
    for( k in 1:20){
      tmp<- timeOffGrid(m,n,np)
      out<- rbind( out, tmp)
    }
  }
  
}

table1<- tapply(out[,3], list( out[,1], out[,2]), mean )
table2<- tapply(out[,4], list( out[,1], out[,2]), mean )

print( xtable( round(table1,3)),file="timingTable1.tex" )
print( xtable( round(table2,4)), file="timingTable2.tex" )

print( xtable( round(cbind(table1,table2),3)),
       file="timingTable3np4.tex" )




np<- 2
out<- NULL
for(  m in c( 128, 256, 512,  2048)){
  for( n in c( 400,  1600, 6400)){
    for( k in 1:20){
      tmp<- timeOffGrid(m,n,np)
      out<- rbind( out, tmp)
    }
  }
  
}

table1<- tapply(out[,3], list( out[,1], out[,2]), mean )
table2<- tapply(out[,4], list( out[,1], out[,2]), mean )


print( xtable( round(cbind(table1,table2),3)),
       file="timingTable3np2.tex" )


np<- 4
out<- NULL
for(  m in c( 128, 256, 512)){
  for( n in c( 400,  1600, 6400)){
      tmp<- timeCEWithNN(m,n,5,np)
      out<- rbind( out, tmp)
      print( tmp)
  }

}
timingCE<-out
save(timingCE,file="timeCE.rda")
out2<- data.frame( 
  cbind( as.integer(out[,1]), as.integer(out[,2]),
            round(out[,4:6],2), 
              
         round(out[,7],3),round(out[,c(8,10)],2) )
)
names( out2)<- c( "m", "n", "CE Setup", "Off Setup","CE", "Off Grid", "predict",
                  "total")
out2<- as.matrix( out2
)
library( xtable)

print( xtable( out2,
               display=c("s","d","d",rep("f",6))
              )
               ,  
 include.rownames=FALSE,
       file="timingTable3np4.tex"
 )
       





