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
circulantEmbeddingSetup <- function( 
     grid, M = NULL, 
     cov.function="stationary.cov", cov.args=NULL,
     delta=NULL, ...) {
    #
    # if cov object is missing then create
    # basically need to enlarge domain and find the FFT of the
    # covariance
    #
    cov.args<-c( cov.args, list(...))
        L<- length( grid)
        dx<- rep( NA, L)
        m<- rep( NA, L)
        for( i in 1:L){
          gridTmp<- grid[[i]]
          dx[i]<- gridTmp[2]- gridTmp[1]
          m[i]<- length( gridTmp)
        }
       
# M is the larger grid  for circulant covariance that includes m 
# should be at least 2*m for embedding to be exact.   
        if( !is.null(delta)){
            M<- rep( NA, L)
            for( i in 1:L){
                M[i]<- m[i] + ceiling( delta/ dx[i])
            }  
        }
# choose a good composite M is not specified        
        if( is.null(M)){
# table of composite M with factors of 2 and 3
            p23<-expand.grid( 0:15, 0:15)
            value<-  2^p23[,1] * 3^p23[,2]
            # all values from 2 and 3 up to 100000
            value<- sort( value[ value <= 1e6])
          M<- rep( NA, L)
          for( i in 1:L){
              # smallest composite choice >= to 2*m
              valueGreater<- value[ value >= 2*m[i] ]
            M[i]<- min( valueGreater )
          }
        }
        if( length(M)!= length( grid)){
          stop("M should be same length as grid")
        }
# create the larger multigrid using M
        bigIndex<- makeMultiIndex( M)
        MCenter<- round( M/2)
        center<-      rbind( MCenter* dx)
# this might be made more efficient another way ....  
        bigGrid<- array( NA, dim(bigIndex) )
        for( i in 1:L){
          bigGrid[,i]<- bigIndex[,i]*dx[i]
        }
        #
        # here is where the actual covariance form is used
        # note passed arguments from call for parameters etc.
        #
        out<- do.call(cov.function, c(cov.args, list(x1 = bigGrid, x2 = center)))  
        # coerce to an array note that this depends on the bigIndex varying in the right way
        out<- array( c(out),M)
        
        #
        # this normalization can be skipped because the simulated field 
        # is stationary and periodic.
        # for example
        # wght <- fft(out)/prod(M)
        # OLD CODE:
        # a simple way to normalize. This could be avoided by
        # translating image from the center ...
        # add to the middle point in the array -- matches the center from above
         temp <- array( 0, M)
         temp[rbind( MCenter)] <- 1
         wght <- fft(out)/(fft(temp) * prod(M))
        
        #
        # wght is the discrete FFT for the covariance suitable for fast
        # multiplication by convolution.
        #
        covObject <- list(m = m, grid = grid, dx=dx, M = M, delta=delta, 
            wght = wght,  call = match.call())
        return( covObject)
       
}
