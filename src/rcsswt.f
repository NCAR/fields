c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
c**********
      subroutine rcsswt(n,y, sy, wt, par)
      implicit double precision (a-h, o-z)
      real*8 y(n), sy(n), wt(n),psi,a,am1,c
      real*8 par(2)
c
c   psi(u) is the derivative of rho(u) defined in rcssr above
c
c   It is composed of peicewise linear and peicewise constant segements
c   and will be continuous except at u for a!=.5. 
c
c
        a= par(2)
        am1 = (1-par(2))
        c= par(1)
        do 100 k=1, n
c   find scaled residual
               r= (y(k)- sy(k))/c
               if( (r.gt. 0)) then 
                         if( r.lt. 1) then 
                               psi= 2*a*r
                               
                         else
                              psi= 2*a
                         endif
               else 
                         if( r.gt.-1) then 
                               psi= 2*am1*r
                               
                         else
                              psi= -2*am1
                         endif

               endif
c
c note weights supplied to cubic spline routine follow the convention that
c they are in terms of standard deviations. ( The more common form is
c   as reciprocal variances
c
        wt(k) = dsqrt( 2*r/psi)
  100   continue
        return
        end
