

addColorBarTriangle <- function(lowerColor=NULL,
                                upperColor=NULL,
                                horizontal=TRUE) {
usr<- par()$usr
plt<- par()$plt

par(xpd=TRUE)
if(horizontal){
  deltaY<-  (usr[4] - usr[3])
  unitX<-   (usr[2] - usr[1])/(plt[2] - plt[1])
  deltaX<-  unitX*(plt[4] - plt[3])
}
else{ 
  deltaX<-  (usr[2] - usr[1])
  unitY<-   (usr[4] - usr[3])/(plt[4] - plt[3])
  deltaY<-  unitY*(plt[2] - plt[1])
}
#
if( !is.null(upperColor) ){
  if(horizontal){
    triangleUpper<- rbind( 
                     c( usr[2],              usr[3] ),
                     c( usr[2]+ deltaX,      usr[3] + deltaY/2 ),
                     c( usr[2],              usr[4])
    )
  }
  else{
    triangleUpper<- rbind( c( usr[1],          usr[4] ),
                       c( usr[1] + deltaX/2,  usr[4] + deltaY ),
                       c( usr[2],          usr[4])
    )
  }
  polygon( triangleUpper, col=upperColor, border=upperColor)
  lines(triangleUpper)
}
#
if(!is.null(lowerColor)){
  if(horizontal){
    triangleLower<- rbind( c( usr[1],          usr[3] ),
                   c( usr[1]- deltaX,  usr[3] + deltaY/2 ),
                   c( usr[1],          usr[4])
    )
  }
  else{
    triangleLower<- rbind( c( usr[1],          usr[3] ),
                           c( usr[1]+ deltaX/2,  usr[3] - deltaY ),
                           c( usr[2],          usr[3])
    ) 
  }
  polygon( triangleLower, col=lowerColor, border=lowerColor)
  lines(triangleLower)
}
#
par(xpd=FALSE)
}
