c fields, Tools for spatial data
c Copyright (C) 2017, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
  
  
  

      subroutine css(h,npoint,x,y,wght,sy,trace,diag,vlam,
     +                 ngrid,xg,yg,job,ideriv,ierr)  
   
C COMPUTES A CUBIC SMOOTHING SPLINE FIT TO A SET OF DATA GIVEN  
C NPOINT(=NUMBER OF DATA VALUES) AND LAMBDA(=VALUE OF  
C THE SMOOTHING parameter). THE CODE IS ADAPTED FROM A  
C PROGRAM IN DEBOOR,C.(1978).A PRACTICAL GUIDE TO SPLINES.  
C SPRINGER-VERLAG : NEW YORK. AN O(NPOINT) ALGORITHM  
C SUGGESTED BY HUTCHINSON AND DEHOOG IS USED TO COMPUTE  
C LVRAGE VALUES AND CONSTRUCT CONFIDENCE INTERVALS.  
c Adapted from Randy Eubank Texas A&M  
c  
c  
c this subroutine solves the problem:  
c  
c   minimize  (1/n) sum(i=1,n) [ (y(i) - f(x(i)))/wght(i) ]**2 + lambda*J(f)  
c    over f  
c   The solution will always be a piecewise cubic polynomial with join   
c   points at the x values. Natural boundary conditions are assumed: at the   
c   x(1) and x(npoints) the second and third derivatives of the slotuion   
c   will be zero   
c   All matrix calculations are done in double precision   
c  
c   Arguments of evss:   
c    h : natural log of lambda ( more convenient scale whepassing a   
c                               real*4)  
c        if h is passed with a value less than or equal -1000 no smoothing will   
c        be done and the spline will interploate the data points  
c    npoint: number of observations  
c    (x,y) : pairs of data points to be smoothed  
c    sy : on return predicted values of f at x  
c    wght : weights used in sum of squares. If the y have unequal   
c      variance then an appropriate choice for wght is the standard deviation  
c      (These weights are not normalized so multiplying by a constant   
c      will imply solving the minimization problem with a different smoothing   
c      parameter)    
c   trace: in matrix from Yhat= A(lambda)Y  with A(lambda) an nxn matrix  
c          trace= tr(A(lambda)) = " effective number of paramters"  
c   diag: diagonal elements of A(lambda) ( this is the mostt   
c         computationally intetnsive opertation in this subroutine.)  
c   vlam: value of the generalized cross-validation functyion (Used to  
c         select an appropriate value for lambda based on the data.)  
c   ngrid,xg,yg : on return, the ith deriv of the spline estimate is  
c                 evaluated on a grid, xg, with ngrid points.  The  
c                 values are returned in yg  
c  
c   ideriv: 0 = evaluate the function according to the job code  
c           1 = evaluate the first derivative of the function according  
c               to the job code  
c           2 = evaluate the second derivative of the function according  
c               to the job code  
c  
c   job: is a vector of three integers
c        (a,b,c)  (a=igcv, b=igrid, c=sorting)
c        a=0  evaluate spline at x values, return predicted values in sy  
c        a=1  same as a=0 plus returning values of trace, diag and vlam  
c        a=2  do no smoothing, interpolate the data  
c        a=3  same as a=1 but use the passed value invlam argument as  
c             a cost an offset factors in the diag vector  
c  
c        b=0  do not evaluate the spline at any grid points  
c        b=1  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)  
c             of the spline at the (unique) sorted,data points, xg, return yg  
c        b=2  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)  
c             of the spline at ngrid equally spaced points between x(1)  
c             and x(npoints)      
c        b=3  evaluate the spline (ideriv=0) or the ith deriv (i=ideriv)  
c             of the spline at ngrid points on grid passed in xg.
c             NOTE: extrapolation of the estimate beyond the range of   
c                     the x's will just be a linear function.  
c  
c
c
c         c=0    X's may not be in sorted order 
c         c=1    Assume that the X's are in sorted order
c         c=2 Do not sort the X's use the  current set of keys
c              Should only be used on a second call to smooth
c              same design
c
c  ierr : on return ierr>0 indicates all is not well :-(
c
      parameter (NMAX=20000)  
      implicit double precision (a-h,o-z)  
      double precision h,trace,vlam  
      double precision wght(npoint),X(npoint),Y(npoint)
      double precision sy(npoint),diag(npoint)  
      double precision xg(ngrid),yg(ngrid)  
      double precision A(NMAX,4),V(NMAX,7)  
      double precision P,SIXP,SIX1MP,cost  
      double precision ux(NMAX),uy(NMAX), uw(NMAX),ud(NMAX),utr 
      integer imx(NMAX) 
      integer idx(NMAX)  
      integer npoint,igcv,igrid, isort,  job(3)  
   
      eps= 1d-7
      cost=1.0  
      offset= 0
        igcv= job(1)
        igrid=job(2)
        isort=job(3)

c  
      if( npoint.gt.NMAX) then  
           ierr=1  
           return  
       endif 
c initialize unique vector and two pointers
      
      do 5 k=1,npoint
      ux(k)=x(k)

      idx(k)=k
    5 continue

c sort vector X along with the keys in imx
c
c       initialize the indices  imx
c
      if( isort.le.1) then 

          do 6 k=1,npoint
            imx(k)=k
    6     continue
       endif
c   
c    sort the X's permuting the indices
c 
       if(isort.eq.0) then 
          call sortm( ux, imx,npoint)
c   the rank of the value x( imx(k)) is k
       endif 
c put y and the weights in the  sorted order 
c  
C**** construct vector consisting of unique observations.  
        
       jj=1 
       ind= imx(1) 
       ux(jj)= x(ind)  
       uy(jj)= y(ind)  
       uw(jj)= wght(ind)  
       idx(1)=jj  
c      normalize eps by the range of the X values
       eps= eps/(  x( imx(npoint)) - x( imx(1)) )
c
c
             do 10 k=2,npoint
c   we are looping through ranks but this is not how the 
c   order of the X are stored. The location of the kth smallest
c   is at idx(k) 
           kshuf= imx(k)  
          if( abs( ux(jj)-x(kshuf)).lt.eps) then  
c**** we have a repeat observation, update weight and the weighted   
c***** average at this point  
            temp1=  1.0d0/uw(jj)**2  
            temp2=  1.0d0/wght(kshuf)**2  
            temp3 = (temp1 + temp2)  
            uy(jj)=  ( uy(jj)*temp1 + y(kshuf)*temp2)/temp3  
            uw(jj)= 1.0d0/dsqrt(temp3)  
          else  
            jj= jj+1  
            ux(jj)= x(kshuf)  
            uy(jj)= y(kshuf)  
            uw(jj)= wght(kshuf)  
          endif
c  save the value that indexes unique values to repliacted ones.
c    x(k) corresponds to the unique X at idx(k)  
          idx(k)=jj  
   10 continue  
      nunq= jj  
  
        itp=0  
        if(igcv.eq.2) itp=1  
  
c   handle condition for interpolation        
  
      if(itp.eq.0) then   
          P=1.d0/(npoint*dexp(h)+ 1.d0)  
      else  
          P=1.d0  
      endif  

      
      call dSETUP(uX,uW,uY,nunq,V,A(1,4),NMAX,itp,ierr)  
C****  check for duplicate X's if so exit   
      if(ierr.gt.0) then  
         return  
      endif  
  
      call dCHOLD(P,V,A(1,4),nunq,A(1,3),A(1,1),NMAX)  
  
c compute predicted values   
  
      SIX1MP=6.d0*(1.d0-P)  
      if(itp.eq.0) then   
         DO 61 I=1,Nunq  
            a(i,1)=uY(I) - SIX1MP*(uW(I)**2)*A(I,1)  
   61    continue  
  
c fill in predicted values taking into account repeated data  
         do 70 k=1,npoint  
            jj= idx(k)  
            sytemp= a(jj,1)
c 
c unscramble the smoothed Y's
           kshuf= imx(k)  

           sy(kshuf)=sytemp  
   70    continue  
      else  
         do 60 i=1,nunq  
            a(i,1)=uy(i)  
   60    continue  
      endif  
  
c return estimates on unique x's if igrid =1  
      if(igrid.eq.1) then  
         do 65 i=1,nunq  
            xg(i)=ux(i)  
             yg(i)=a(i,1)  
   65    continue  
         ngrid=nunq  
      endif  
      if(igrid.ge.1) then   
c  
c********* evaluate spline on grid  
C piecewise cubic COEFFICIENTS ARE STORED IN A(.,2-4).   

         SIXP=6.d0*P  
         DO 62 I=1,nunq  
          A(I,3)=A(I,3)*SIXP
  62     continue
         NPM1=nunq - 1  
         DO 63 I=1,NPM1  
            A(I,4)=(A(I+1,3)-A(I,3))/V(I,4)  
            A(I,2)=(A(I+1,1)-A(I,1))/V(I,4)  
     *        -(A(I,3)+A(I,4)/3.*V(I,4))/2.*V(I,4)  
   63    continue  
c  
c  create equally spaced x's if asked for ( igrid=2)  
c  
         if( igrid.eq.2) then   
              step= (ux(nunq)-ux(1))/(ngrid-1)  
              xg(1)=ux(1)  
              do 190 j=2,ngrid-1  
                 xg(j)= xg(j-1)+step  
  190         continue  
              xg(ngrid)=ux(nunq)  
         endif  
  
  
         uxmin= ux(1)  
         uxmax= ux(nunq)  
  
         a1= a(1,1)  
         an= a(nunq,1)  
  
         b1= a(1,2)  
  
         d= ux(nunq)- ux(nunq-1)  
         ind= nunq-1  
         bn= a(ind,2) + a(ind,3)*d + a(ind,4)*(d**2)/2.d0  
  
c evalute spline by finding the interval containing the evaluation   
c point and applying the cubic formula for the curve estiamte  
c finding the interval such that ux(ind) <=xg(j) < ux( ind+1)  
c is done using a bisection search   
  
         do 195 j=1,ngrid  
           lint= ifind(xg(j),ux,nunq)  
          if( lint.eq.0) then   
          d= xg(j)-uxmin  
               
            if (ideriv .eq. 0)   
     -         yg(j)= a1 + b1*d   
            if (ideriv .eq. 1)  
     -         yg(j)=  b1  
            if (ideriv .eq. 2)  
     -         yg(j)= 0.0  
         endif  
          if( lint.eq.nunq) then   
          d= xg(j)-uxmax  
               
            if (ideriv .eq. 0)   
     -         yg(j)= an + bn*d   
            if (ideriv .eq. 1)  
     -         yg(j)=  bn  
            if (ideriv .eq. 2)  
     -         yg(j)= 0.0  
         endif  
            if( ((lint.ge.1 ) .and.( lint.lt.nunq))) then   
                ind=lint  
c            a1=a(ind,1)  
c            a2=a(ind,2)  
c            b=a(ind,3)/2.d0  
c            c=a(ind,4)/6.d0  
c  
            d= xg(j)-ux(ind)  
  
            if (ideriv .eq. 0)   
     -         yg(j)= a(ind,1) + a(ind,2)*d + a(ind,3)*(d**2)/2.d0  
     -                + a(ind,4)*(d**3)/6.d0  
            if (ideriv .eq. 1)  
     -         yg(j)= a(ind,2) + a(ind,3)*d + a(ind,4)*(d**2)/2.d0  
            if (ideriv .eq. 2)  
     -         yg(j)= a(ind,3) + a(ind,4)*d  
         endif  
c  
  195  continue   
      endif  
c****end of evaluation block  
  
      if((igcv.eq.1).or.( igcv.eq.3)) then  
        if( igcv.eq.3)  then 
              cost= diag(1)
              offset=diag(2)
         endif
c*****                       begin computing gcv and trace  
C COMPUTE LVRAGE VALUES ,THE VARIANCE ESTIMATE   
C     SGHAT2 AND CONFIDENCE INTERVALS.  
c  
         call dLV(nunq,V,uw,SIX1MP,utr,ud,NMAX)  
  
         rss=0.d0  
         trace=0.d0  
         vlam=0.d0  
  
         do 100 i=1,nunq  
c            rss= rss + ((uy(i)-a(i,1))/uw(i))**2  
            trace= trace +ud(i)  
  100    continue  
  
         do 110 k=1,npoint 
            kshuf= imx(k)  
            jj= idx(k) 
            diag(kshuf)= ud(jj)
            rss= rss + ( (y(kshuf)- sy(kshuf))/wght(kshuf) )**2 
  110    continue  
         ctrace=  2+  cost*( trace-2)   
         if( (npoint -ctrace -offset) .gt. 0.d0) then  
            vlam= (rss/npoint)/( 1- (ctrace-offset)/npoint)**2   
         else  
            vlam=1e20  
         endif  
  
      endif  
      
       
  
      return  
  
      END  
