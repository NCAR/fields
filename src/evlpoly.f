c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 


      subroutine evlpoly(x,n,coef,j,result)

c evaluates a polynomial: coef(1) + sum_i= 2,j coef(i)x**i

      integer j,n, i
      double precision x(n), result(n), coef(j)
      double precision temp, tempx, temp2
      
      do 10 i = 1,n

         temp= coef(1)
         tempx= x(i)
         temp2= tempx 

c      temp is set to constant now loop over powers

         do 20 kk= 2, j
           temp= temp + coef(kk)* temp2
           temp2= temp2*tempx
   20    continue

      result(i)= temp

   10 continue
 
      return
      end
