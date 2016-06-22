c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 

c
c evaluates exponential radial bais function
c 
       subroutine multeb( nd,x1,n1, x2,n2, par, c,h,work)
       implicit double precision (a-h,o-z)
       integer nd,n1,n2,ic
       
       real*8 par(nd),x1(n1,nd), x2(n2,nd), c(n2), h(n1),sum
       real*8 work( n2), ddot

c****** work aray must be dimensioned to size n2
c **** loop through columns of output matrix K
c*** outer most loop over columns of x1 and x2 should reduce paging 

       do 5 ir= 1, n1

c
 
c evaluate all basis functions at  x1(j,.)       
       do 10 j =1,n2
c
c  zero out sum accumulator
c
         sum=0.0
      do 15  ic=1,nd
c
c** accumulate squared differences
c 

            sum= sum+ (dabs(x1(ir,ic)- x2(j,ic)))**2

 15             continue
        work(j)=sum
 10    continue

C**** evaluate squared distances  with basis functions. 

          call expfn( n2,work(1),par)
c
c***** now the dot product you have all been waiting for!
c
          h(ir)= ddot( n2, work(1), 1, c(1),1)
 5      continue

       return
       end
