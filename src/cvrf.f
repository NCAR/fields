c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 


c** finds value of h minimizing the  generalized cross-validation

      subroutine cvrf
     * (h,nobs,x,y,wt,sy,trace,diag,din,dout,cv,ierr)  
      implicit double precision (a-h,o-z)    
      REAL*8 h,trace,cv
      REAL*8 x(nobs),y(nobs),wt(nobs)
      REAL*8 sy(nobs),diag(nobs),dumm1(1),dumm2(1)
      real*8 din(10), dout(10)
      integer ngrid, ierr, job(3),ideriv
      data job/3,0,0/
      data ideriv,ngrid/0,1/
      call rcss(h,nobs,x,y,wt,sy,trace,diag,cv,  
     +          ngrid,dumm1,dumm2,job,ideriv,din,dout,ierr)  
      nit= int( dout(1))
      trace=dout(3)
           
c
      return
      end
