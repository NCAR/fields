c fields, Tools for spatial data
c Copyright (C) 2017, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
      subroutine multWendlandG( mx, my, deltaX, deltaY,
     *  nc, center, coef, h, flag)
      integer mx, my, nc, flag
      double precision  deltaX, deltaY, center(nc,2), coef(nc)
      double precision h(mx,my) 
c 
      integer j, k, l, m1, m2, n1, n2 
      double precision  kstar, lstar, d, wendlandFunction
      do j = 1, nc
         kStar= center(j,1)
         lStar= center(j,2)
         m1 =  max( ceiling(-deltaX + kStar),  1)
         m2 =  min(   floor( deltaX + kStar), mx)
         n1 =  max( ceiling(-deltaY + lStar),  1)
         n2 =  min(   floor( deltaY + lStar), my)
         do l = n1, n2
            do k = m1, m2
              d = dsqrt( ((k-kStar)/deltaX)**2 + ((l- lStar)/deltaY)**2)
              h(k,l) = h(k,l) + wendlandFunction( d)* coef(j)
c             h(k,l) = h(k,l) + 2
            enddo
         enddo
      enddo
      flag=0
      return
      end

      double precision function wendlandFunction(d)
      double precision d
         if( d.GE.1) then 
           wendlandFunction = 0
         else
           wendlandFunction = ((1-d)**6) * (35*d**2 + 18*d + 3)/3
         endif
      return
      end
