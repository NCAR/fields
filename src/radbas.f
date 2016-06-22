c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
c**** subroutine to fill in the omega ( or K) matrix for 
c**** ridge regression S funcion
c**** K_ij= radfun( distance( x1_i, x2_j))
c
       subroutine radbas( nd,x1,n1, x2,n2, par, k)
       integer nd,n1,n2,ic   
       double precision par(2), x1(n1,nd), x2(n2,nd), k(n1,n2)
       double precision xtemp, radfun
c **** loop through columns of output matrix K
c*** outer most loop over columns of x1 and x2 should reduce memory swaps 
       do  ic= 1, nd
          do  j =1,n2
              xtemp= x2(j,ic)
              do   i= 1, n1
c** accumulate squared differences
                 k(i,j)=  (x1(i,ic)- xtemp)**2 + k(i,j)
              enddo  
          enddo
       enddo
c**** at this point k( i,j) is the squared distance between x1_i and x2_j
c*** now evaluate radial basis functions
       do  j =1,n2
          do  i= 1, n1
              k(i,j)=  radfun( k(i,j), par(1), par(2) )
          enddo
       enddo
       return
       end
