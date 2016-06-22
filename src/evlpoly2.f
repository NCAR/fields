c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 


      subroutine evlpoly2(x,n,nd,ptab,j,coef,result)

c evaluates a polynomial: coef(1) + sum_i= 2,j coef(i)x**i

      integer j,n, i, nd
      double precision x(n,nd), result(n), coef(j)
      integer ptab(j,nd)
      double precision temp,  temp2
      
      do 10 i = 1,n

         temp= 0 

c for a given vector accumlate the polynomial terms

         do 20 kk= 1, j
            temp2 =1.0

            do 30 l=1, nd
                   if( ptab(kk,l).ne.0) then
                   temp2= temp2* (x(i,l)**ptab(kk,l))
                   endif
 30         continue

          temp = temp + temp2*coef(kk) 

   20    continue

      result(i)= temp

   10 continue
 
      return
      end
