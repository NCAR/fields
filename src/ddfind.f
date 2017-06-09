c fields, Tools for spatial data
c Copyright (C) 2017, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
c
       subroutine ddfind( nd,x1,n1, x2,n2, D0,ind,rd,Nmax, iflag)

       integer nd,n1,n2, ind(nmax,2)
       integer kk, i,j, ic
       
       double precision x1(n1,nd), x2(n2,nd), D0, rd(Nmax), D02, dtemp

c****   counter  for accumulating close points
        kk=0 
        D02= D0**2
            do  15 i= 1, n1
                  do 10 j =1,n2

c
c** accumulate squared differences
c
               dtemp= 0.0
               do 5 ic= 1, nd 
                    dtemp= dtemp + (x1(i,ic) - x2(j,ic))**2
                    if( dtemp.gt.D02) goto 10
 5             continue
               
c****       dtemp is less than D0 so save it as a close point

             kk=kk+1

c**** check if there is still space 
             if( kk .gt. Nmax) then 
                iflag= -1
                goto 20 
             else
                ind(kk,1)= i
                ind(kk,2)= j
                rd(kk)= sqrt( dtemp)
             endif
       

 10             continue
 15         continue
 
         Nmax=kk  
20       continue


       return
       end
