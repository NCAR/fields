c fields, Tools for spatial data
c Copyright (C) 2017, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
       

C** evaluates partial derivatives of radial basis functions with
c** nodes at x2 and at the points x1
c
       subroutine mltdrb( nd,x1,n1, x2, n2, par, c,h,work)
       implicit double precision (a-h,o-z)
       integer nd,n1,n2,ic, ivar, ir, j     
       double precision par(2), x1(n1,nd)
       double precision  x2(n2,nd), c(n2), h(n1, nd)
       double precision sum, sum2
       double precision work( n2)
c      double precision ddot
       do 1000 ivar=1, nd
c****** work aray must be dimensioned to size n2
c **** loop through columns of output matrix K
c***  outermost loop over columns of x1 and x2 should
c***  help to access adjacent values in memory.          
       do 5 ir= 1, n1
c evaluate all basis functions at  x1(j,.)       
       do 10 j =1,n2
c  zero out sum accumulator
          sum=0.0
          do 15  ic=1,nd
c** accumulate squared differences
             sum= sum+ (x1(ir,ic)- x2(j,ic))**2
 15       continue
          work(j)=sum
 10    continue
C**** evaluate squared distances  with basis functions. 
       call drdfun( n2,work(1) ,par(1) )
       do 11 j= 1, n2
          work( j)= 2.0*work(j)*(x1(ir,ivar)- x2(j,ivar))
 11    continue
c*****now the dot product you have all been waiting for!
       sum2= 0.0
       do 12 j = 1, n2
          sum2 = sum2  + work(j)*c(j)
 12    continue
c     h(ir,ivar)= ddot( n2, work(1), 1, c(1),1)
       h(ir,ivar) = sum2
 5     continue
 1000   continue
       return
       end
