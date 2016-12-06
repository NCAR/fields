c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
  
      subroutine csstr(h,nobs,x,y,wght,c,offset,trace,vlam,work,ierr)
      parameter(mxM=20000)
      implicit double precision (a-h,o-z)
      double precision h,trace, vlam,c,offset
      double precision x(nobs),y(nobs),wght(nobs)
      double precision work(nobs),diag(mxM),dumm1(1),dumm2(1)
      integer job(3),ideriv,ierr, ndum
      data ideriv/0/
       job(1)=3
       job(2)=0
       job(3)=0
       diag(1)=c
       diag(2)=offset
       ndum=1
      call css(h,nobs,x,y,wght,work,trace,diag,vlam,ndum,dumm1,dumm2,
     -            job,ideriv,ierr)

      return
      end
