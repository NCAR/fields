c fields, Tools for spatial data
c Copyright (C) 2017, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 


       subroutine mltdtd( nd,x1,n1,np,ptab, d,h)
       implicit double precision (a-h,o-z)
       integer nd,n1,np, ivar
       
       double precision x1(n1,nd)
       double precision   d(np), h(n1, nd)
       double precision work,  prod, xp,xx
       integer ptab(np,nd)
c  outer most loop is over the variables w/r partial derivative

       do 1000 ivar=1, nd
c****** work aray must be dimensioned to size np
       
       do 5 ir= 1, n1
c next loop is over rows of x1

c
 
c evaluate all partials of polynomials  at  x1(j,.)       
c take ddot product of this vector with d 
c this is the element to return in h(ir,ivar)
       work=0.0
       do 10 j =1,np 
        prod=0.0
       ipv= ptab( j,ivar)
        if( ipv.gt.0) then
            prod=1.0    
            do 11 k= 1, nd
               ip= ptab(j,k)
c ip is the power of the kth variable in the jth poly
             if( ip.eq.0) goto 11 
                 xx= x1(ir,k)
c**
               if( k.eq.ivar) then
                   if( ip.eq.1) then 
                        xp=1.0
                   else 
                        xp= (ip)* xx**(ip-1) 
                  endif
               else 
                  xp= xx**(ip)
               endif
                prod=prod*xp
 11            continue
        endif
        work= work + prod* d(j)
10      continue

c
        h(ir,ivar)=work
 5      continue
 1000   continue
       return
       end
