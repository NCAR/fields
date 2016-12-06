c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
      SUBROUTINE dSETUP(X,WGHT,Y,NPOINT,V,QTY,NMAX,itp,ierr)
C   PUT DELX=X(.+1)-X(.) INTO V(.,4)
C   PUT THE THREE BANDS OF THE MATRIX Q-TRANSP*D INTO
C     V(.,1-3)
C   PUT THE THREE BANDS OF (D*Q)-TRANSP*(D*Q) AT AND 
C     ABOVE THE DIAGONAL INTO V(.,5-7)
C   HERE Q IS THE TRIDIAGONAL MATRIX OF ORDER (NPOINT
C     -2,NPOINT) THAT SATISFIES Q-TRANSP*T=0 AND WGHT
C     IS THE DIAGONAL MATRIX WHOSE DIAGONAL ENTRIES 
C     ARE THE SQUARE ROOTS OF THE WEIGHTS USED IN THE 
C     PENALIZED LEAST-SQUARES CRITERION
c
      implicit double precision (a-h,o-z)
      double precision WGHT(NMAX),X(NMAX),y(NMAX)
      double precision QTY(NMAX),V(NMAX,7)
      double precision DIFF,PREV
      INTEGER NPOINT,I,NPM1
c
      NPM1=NPOINT -1
      V(1,4)=X(2)-X(1)
      if(v(1,4).eq.0.d0) then
                ierr=5
                return
      endif

      DO 11 I=2,NPM1
         V(I,4)=X(I+1) - X(I)
         if(v(I,4).eq.0.d0) then
                ierr=5
                return
         endif
         if(itp.eq.0) then 
            V(I,1)=WGHT(I-1)/V(I-1,4)
            V(I,2)=-WGHT(I)/V(I,4) - WGHT(I)/V(I-1,4)
            V(I,3)=WGHT(I+1)/V(I,4)
         else
            V(I,1)=1.d0/V(I-1,4)
            V(I,2)=-1.d0/V(I,4) - 1.0/V(I-1,4)
            V(I,3)=1.d0/V(I,4)
         endif
   11 continue
c
      V(NPOINT,1)=0.d0
      DO 12 I=2,NPM1
        V(I,5)=V(I,1)**2 + V(I,2)**2 + V(I,3)**2
 12   continue
      IF(NPM1 .LT. 3)GO TO 14
      DO 13 I=3,NPM1
         V(I-1,6)=V(I-1,2)*V(I,1) + V(I-1,3)*V(I,2)
 13      continue
   14 V(NPM1,6)=0.d0
      IF(NPM1 .LT. 4)GO TO 16
      DO 15 I=4,NPM1
       V(I-2,7)=V(I-2,3)*V(I,1)
 15   continue   
   16 V(NPM1-1,7)=0.d0
      V(NPM1,7)=0.d0
c
C  CONSTRUCT Q-TRANSP. * Y IN QTY
      PREV=(Y(2) - Y(1))/V(1,4)
      DO 21 I=2,NPM1
         DIFF=(Y(I+1)-Y(I))/V(I,4)
         QTY(I)=DIFF - PREV
         PREV=DIFF
 21   continue
c
      RETURN
      END
