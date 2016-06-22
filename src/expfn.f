c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 

       subroutine expfn(n,d2, par)
       real*8 d2(n), par(1)
       integer n

         do 5 k =1,n
         d2(k)=  exp(-1*d2(k)**(par( 1)/2))
   5     continue

        return
        end
