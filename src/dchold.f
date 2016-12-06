c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 




      SUBROUTINE dCHOLD(P,V,QTY,NPOINT,U,QU,NMAX)
c CONSTRUCTS THE UPPER THREE DIAGONALS IN V(I,J),I=2,
C   NPOINT-1,J=1,3 OF THE MATRIX 6*(1-P)*Q-TRANSP*
C   (D**2)*Q + P*R . (HERE R IS THE MATRIX PROPORTIONAL
C   TO Q-TRANSP*KSUBN*Q , WHERE KSUBN IS THE MATRIX 
C   WITH ELEMENTS K(X(I),X(J)) AND K IS THE USUAL
C   KERNEL ASSOCIATED WITH CUBIC SPLINE SMOOTHING.)
C   THE CHOLESKY DECOMPOSITION OF THIS MATRIX IS COMPUTED
C   AND STORED IN V(.,1-3) AND THE EQUATION SYSTEM FOR 
C   THE QUADRATIC COEFFICIENTS OF THE SPLINE IN ITS 
C   PIECEWISE POLYNOMIAL REPRESENTATION IS SOLVED . THE 
C   SOLUTION IS RETURNED IN U.
c
      double precision P,QTY(NMAX),QU(NMAX),U(NMAX),V(NMAX,7)
      double precision SIX1MP,TWOP,RATIO,PREV
      INTEGER NPOINT,I,NPM1,NPM2
c
      NPM1=NPOINT - 1
C****  CONSTRUCT 6*(1-P)*Q-TRANSP.*(D**2)*Q + P*R
      SIX1MP=6.d0*(1.d0 - P)
      TWOP=2.d0*P
      DO 2 I=2,NPM1
         V(I,1)=SIX1MP*V(I,5) + TWOP*(V(I-1,4)+V(I,4))
         V(I,2)=SIX1MP*V(I,6) + P*V(I,4)
         V(I,3)=SIX1MP*V(I,7)
   2  continue
      NPM2=NPOINT - 2
      IF(NPM2 .GE. 2)GO TO 10
      U(1)=0.d0
      U(2)=QTY(2)/V(2,1)
      U(3)=0.d0
      GO TO 41
C  FACTORIZATION : FACTORIZE THE MATRIX AS L*B-INV*
C     L-TRANSP , WHERE B IS A DIAGONAL MATRIX AND L
C     IS UPPER TRIANGULAR.
   10 DO 20 I=2,NPM2
         RATIO=V(I,2)/V(I,1)
         V(I+1,1)=V(I+1,1) - RATIO*V(I,2)
         V(I+1,2)=V(I+1,2) - RATIO*V(I,3)
         V(I,2)=RATIO
         RATIO=V(I,3)/V(I,1)
         V(I+2,1)=V(I+2,1) - RATIO*V(I,3)
         V(I,3)=RATIO
 20   continue
C  FORWARD SUBSTITUTION
      U(1)=0.d0
      V(1,3)=0.d0
      U(2)=QTY(2)
      DO 30 I=2,NPM2
         U(I+1)=QTY(I+1) - V(I,2)*U(I) -V(I-1,3)*U(I-1)
 30   continue
C  BACK SUBSTITUTION
      U(NPOINT)=0.d0
      U(NPM1)=U(NPM1)/V(NPM1,1)
      DO 40 I=NPM2,2,-1
         U(I)=U(I)/V(I,1) - U(I+1)*V(I,2) - U(I+2)*V(I,3)
 40   continue
C  CONSTRUCT Q*U
   41 PREV=0.d0
      DO 50 I=2,NPOINT
         QU(I)=(U(I)-U(I-1))/V(I-1,4)
         QU(I-1)=QU(I) - PREV
         PREV=QU(I)
 50   continue
      QU(NPOINT)=-QU(NPOINT)
c
      RETURN
      END
