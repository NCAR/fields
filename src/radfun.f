c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
c evaluates thin plate spline radial basis function
       double precision function radfun(d2, par1, par2)
       double precision d2, par1, par2 
       if( d2.lt.1e-20) then 
           d2= 1e-20
       endif
       if( int(par2).eq.0) then
           radfun= (d2)**( par1)
       else 
c note: d2 is squared distance
c divide by 2 to have log evaluated on distance
c as opposed to squared distance. 
           radfun= (log(d2)/2) * ((d2)**( par1))
       endif
       return
       end

