c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 






       subroutine rcss(h,npoint,x,y,wt,sy,trace,diag,cv,  
     +                  ngrid,xg,yg,job,ideriv,din,dout,ierr)  
c This a program to compute a robust univariate spline according to the
c model:
c   minimize   (1/n) sum(i=1,n)[rho( y(i)-f(x(i) )] + lambda*J(f)
c    over f
c definition of the rho function and its derivative are in rcsswt
c and rcssr
c
c One way of speeding convergence is to use the results from a previous
c estimate as the starting values for the another estimate. This is
c particulary appropriate when lambda or a parameter in the rho function
c is being varied. Moderate changes in lambda will often yield similar
c estimates. The way to take advantage of this is to pass the weights
c from the previous fit as teh starting values for the next estimate
c   
c   Arguments of rcss:   
c    h : natural log of lambda  
c                               
c        if h is passed with a value less than or equal -1000 no smoothing will   
c        be done and the spline will interploate the data points  
c    npoint: number of observations  
c    (x,y) : pairs of data points to be smoothed  
c            x(i) are assumed to be increasing. Repeated   
c            observations at the same x are not handled by this routine.  
c    sy : on return predicted values of f at x  
c    wt : weights used in the iterivatively reweighted least 
c           squares algorithm to compute robust spline. Vector
c           passed are used as the starting values. Subsequent
c           iterations compute the weights by a call to  the 
c           subroutine     rcsswt  
c
c          that is the linear approximation of teh estimator at 
c              convergence. 
c          trace= tr(A(lambda)) = " effective number of paramters"  
c   diag: diagonal elements of A(lambda) ( this is the most   
c         computationally intetnsive opertation in this subroutine.)  
c   cv: approximate cross-validation function 
c         using the linear approximation at convergence
c 
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
c   din:  Vector of input parameters 
c           din(1)= cost for cv
c           din(2)= offset for cv
c           din(3)= max number of iterations
c           din(4)= tolerance criterion for convergence
c
c           din(5)= C scale parameter in robust function (transition from
c                   quadratic to linear.
c           din(6)= alpha   1/2 slope of rho function for large, positive
c                   residuals   slope for residuals <0 : is  1/2 (1-alpha) 
c                  see comments in rcsswt for defintion of the  rho and psi 
c                  functions
c
c   job: in decimal job=(a,b,c)  (a=igcv, b=igrid)   
c        a=0  evaluate spline at x values, return predicted values in sy  
c        a=1  same as a=0 plus returning values of trace, diag and cv  
c        a=2  do no smoothing, interpolate the data  
c        a=3  same as a=1 but use the passed values in din(1) din(2)
c             for computing cv function
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
c         c=1    Assume that the X's are in sorted order
c         c=2 Do not sort the X's use the  current set of keys
c              Should only be used on a second call to smooth
c              same design
 
c Arguments of subroutine:
c    dout(1)= numit
c    dout(2)=tstop
c    dout(3) = trace
c    dout(4)= cv
c      numit: actual number of iterations for convergence
c      tstop: value of convergence criterion at finish.
c
c ierr: if ierr>0 something bad has happened
c       ierr>10 indicates problems in the call to the cubic spline
c       routine.
c  
   
      parameter(NMAX=20000)  
      implicit double precision (a-h,o-z)  
      REAL*8 h,trace,cv  
      REAL*8 wt(npoint),X(npoint),Y(npoint),sy(npoint),diag(npoint)  
      REAL*8 xg(ngrid),yg(ngrid)  
      REAL*8 din(10), dout(10),cost,  offset, dum1, dum2
      
      integer npoint,ngrid ,itj(3), job(3)


      if(npoint.gt.NMAX) then
           ierr=1
           return
      endif

       maxit= int(din(3))
       tstop=din(4)
       ybar=0.0
       ysd=0.0       
       do 10 k=1,npoint
            diag(k) = y(k)
            ybar= ybar+ y(k)
            ysd= ysd+ y(k)**2
   10  continue 
        ybar= ybar/npoint
        ysd= sqrt( ysd/npoint - ybar**2)     
c Start iterating
          test=0.0

       itj(1)= 0
       itj(2)=0
       itj(3)=0
      do 500 it=1,maxit
      if( it.gt.1) then
        itj(3)=2
      endif
c fit a weighted least squares cubic smoothing spline
      call css(h,npoint,x,y,wt,sy,
     *  dum1,diag,dum2,ngrid,xg,yg,
     *  itj,ideriv,ierr)
c check for an error returned by spline routine
      if(ierr.gt.0) then
c         add 10 so these can be distinguished from errors in this routine
          ierr= 10 + ierr
          return
       endif

c compute convergence criterion
c The intent of this criterion is to find out when successive iterations
c produce changes that are small 


          do 25 i=1,npoint
               test=(diag(i)-sy(i))**2 + test 
               diag(i)= sy(i)
 25       continue

           test=sqrt(test/npoint)/ysd
           if( test.lt.tstop  ) then
c            * exit loop *
               numit=it
               goto 1000
            endif
             
c
c    make up new set of weights
c
          call rcsswt( npoint, y,sy,wt, din(5))
 480       format( 5e15.6)
c    reinitialize test criterion for convergence
c
      test=0.0
500   continue

      numit= maxit
 1000 continue
c One last call if job code is not 0 
      if( (job(1).ne.0).or.(job(2).ne.0)) then
c
      call css(h,npoint,x,y,wt,sy,
     *  trace,diag,cv,ngrid,xg,yg,
     *  job,ideriv,ierr)

      ia= job(1)
      ig= job(2)
c      if(ig.gt.0)  then
c      endif
c calculate cv value if asked for 
      if( (ia.eq.1) .or.( ia.eq.3) ) then 
             if(ia.eq.3) then 
             cost= din(1)
             offset= din(2)/npoint
             else
             cost=1
             offset= 0
             endif

             cv=0.0
             do 1500 k=1,npoint
c compute approx. cross-validated residual
             

c plug cv residual into rho function, din(5) is the begining of parameter 
c vector for the rho function (scale and alpha)
c
c but only do this if the leverage is away from one. 
c this prevents the numerical problems with quantile splein of a zero
c  residual and 
c a zero denominator. 
                 if(1- diag(k).gt.1e-7) then  
                  resid= (y(k)- sy(k))/( 1- cost*(diag(k)+offset))
                  cv= cv + rcssr(resid, din(5))
      endif
    
 1500        continue
             cv= cv/npoint
      endif 
      endif

      dout(1)=numit
      dout(2)=test
      dout(3)=trace
      dout(4)=cv

      return
      end
