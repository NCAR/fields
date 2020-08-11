# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2018
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
"circulantEmbedding" <- function(obj) {
  if( is.null(obj$mDim )){
# backwards compatibility for 2D case     
    mDim<-c( obj$m, obj$n)
    MDim<-c( obj$M, obj$N)
  }
  else{
    mDim<- obj$mDim
    MDim<- obj$MDim
  }
  
    nZ<- length(MDim)
    if (any(Re(obj$wght) < 0)) {
        stop("FFT of covariance has negative\nvalues")
    }
    prodM<- prod( MDim)
    Z <- fft( array(rnorm(prodM), MDim) )
    # create the multindex that is 1:mDim[1] , 1:nDim[2] ..., 1:mDim[,nZ ]
    # two special cases, 1D and 2D  added for efficiency and clarity ...
     if( nZ ==1){ 
       bigIndex<- cbind( 1 : mDim[1] )
     }
    if( nZ==2){
      bigIndex<-  cbind( rep( 1:mDim[1], mDim[2]),  
                         rep( 1:mDim[2], rep(mDim[1], mDim[2]) )
                       )
    }
    if( nZ > 2){
      bigIndex<-  cbind( rep( 1:mDim[1], mDim[2]),  rep( 1:mDim[2], mDim[1]) )
    for ( k in 3 : nZ){
      bigN<- nrow( bigIndex)
      bigIndex<- cbind(
        rep( bigIndex,mDim[k]), 
        rep( 1 : mDim[k], rep( bigN,mDim[k]) )
      )
    }
    }
    
    #bigIndex<- as.matrix(expand.grid( 1:mDim[1], 1:mDim[2]))
    out<-  Re(
               fft(sqrt(obj$wght) * Z, inverse = TRUE)
             )/sqrt(prodM)
    out<- array( out[bigIndex], mDim)
    return( out)
}
