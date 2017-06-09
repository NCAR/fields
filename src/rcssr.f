c fields, Tools for spatial data
c Copyright (C) 2017, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
 
      double precision function rcssr(r,par)
c
c     robust rho function:
c  This is a peicewise polynomial with knots at -C , 0 and C
c  the function is quadratic for -C<u<0 and 0<u<C
c  the function is linear for u<-C and u>C
c  rho is continuous for all u aqnd differentiable for all points 
c  except u=0 when a != 1/2 
c   
c
c    rho(u) =      2*a*u - a*c     for u>C
c                  a*u**2/C        for   0<u< C   
c                  (1-a)*u**2/C    for -C<u<0
c                  2*(1-a)*u - (1-a)*C  for u< -C
c
c        Note a= par(1), C= par(2)
      implicit double precision (a-h, o-z)
      double precision r, par(2),c,a
      c= par(1)     
      if( r.gt.0 ) then 
         a=par(2)
       else
         a =(1-par(2))
         r= -r
      endif 
      if( r.le.c) then
            rcssr= a*r*r/c
      else
           rcssr= 2*a*(r) - a*c
      endif
      return
      end
