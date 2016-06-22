c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 

       subroutine radfn2(n,d2, par)
c derivative of  radial basis function
c  r(d2)=  d2**(par(1)   or  par(2).ne. 0   d2**par(1) log( d2)  
c
       double precision d2(n), par(2), dtemp, p1,p2
       integer n
       if( int(par(2)).eq.0) then

         p1= par(1)   
         p2= par(1)-1   
c NOTE p2 should be nonnegative   
         do  k =1,n
           dtemp= d2(k)
           if( dtemp.lt.1e-20) then
            d2(k)=0.0
           else 
            d2(k)= p1*(dtemp)**( p2)
           endif
         enddo
       else 
         do  k=1,n
          dtemp= d2(k)
          if( dtemp.gt.1e-20)  then
           d2(k)=  p1*log(dtemp)*(dtemp)**(p2) +
     *                dtemp**(p2)
          else
           d2(k)=0.0
          endif
         enddo 
       endif
        return
        end

