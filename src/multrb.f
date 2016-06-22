c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
C** evaluates radial basis functions 
c**** K_ij= radfun( distance( x1_i, x2_j))
c**** and does the multplication  h= Kc
c****  K is n1Xn2
c****  h is n1Xn3
c***** c is n2xn3
 

       subroutine multrb( nd,x1,n1, x2,n2, par, c,n3,h,work)
       implicit double precision (a-h,o-z)
       integer nd,n1,n2,n3,ic, jc,j   
       double precision par(2),x1(n1,nd), x2(n2,nd), c(n2,n3), h(n1,n3) 
       double precision work( n2), ddot, sum
       double precision radfun

c****** work aray must be dimensioned to size n1
c **** loop through columns of output matrix K
      do 5 ir= 1, n1
          do 10 j =1,n2
             sum=0.0
             do 15  ic=1,nd
c** accumulate squared differences
                sum= sum+ (x1(ir,ic)- x2(j,ic))**2
 15          continue
             work(j)=radfun( sum, par(1), par(2))
 10       continue
c***** dot product for matrix multiplication
          do 30 jc=1,n3
             h(ir,jc)= ddot( n2, work, 1, c(1,jc),1)
30        continue  
 5    continue
      return
      end
