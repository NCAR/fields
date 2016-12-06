c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 




      SUBROUTINE dLV(NPOINT,V,WGHT,SIX1MP,TR,LEV,NMAX)
c CONSTRUCTS THE UPPER THREE DIAGONALS OF (6*(1-P)*
C   Q-TRANSP*(D**2)*Q + P*R)-INV USING THE RECURSION
C   FORMULA IN HUTCHINSON,M.F. AND DEHOOG,F.R.(1985).
C   NUMER. MATH. 47,99-106, AND STORES THEM IN V(.,5-7).
C   THESE ARE USED IN POST AND PRE MULTIPLICATION BY
C   Q-TRANSP AND Q TO OBTAIN THE DIAGONAL ELEMENTS OF
C   THE HAT MATRIX WHICH ARE STORED IN THE VECTOR LEV.
C     THE TRACE OF THE HAT MATRIX IS RETURNED IN TR.
c
      double precision V(NMAX,7),TR,W1,W2,W3,SIX1MP
      double precision wght(NMAX)
      double precision LEV(npoint)
      INTEGER NPM1,NPM2,NPM3,NPOINT
c
      NPM1=NPOINT - 1
      NPM2=NPOINT - 2
      NPM3=NPOINT - 3
C   RECURSION FOR DIAGONALS OF INVERSE MATRIX
      V(NPM1,5)=1/V(NPM1,1)
      V(NPM2,6)=-V(NPM2,2)*V(NPM1,5)
      V(NPM2,5)=(1/V(NPM2,1)) - V(NPM2,6)*V(NPM2,2)
      DO 10 I=NPM3,2,-1
         V(I,7)=-V(I,2)*V(I+1,6) - V(I,3)*V(I+2,5)
         V(I,6)=-V(I,2)*V(I+1,5) - V(I,3)*V(I+1,6)
         V(I,5)=(1/V(I,1))- V(I,2)*V(I,6) - V(I,3)*V(I,7)
   10 CONTINUE
C   POSTMULTIPLY BY (D**2)*Q-TRANSP AND PREMULTIPLY BY Q TO
C     OBTAIN DIAGONALS OF MATRIX PROPORTIONAL TO THE 
C     IDENTITY MINUS THE HAT MATRIX.
      W1=1.d0/V(1,4)
      W2= -1.d0/V(2,4) - 1.d0/V(1,4)
      W3=1.d0/V(2,4)
      V(1,1)=V(2,5)*W1
      V(2,1)=W2*V(2,5) + W3*V(2,6)
      V(2,2)=W2*V(2,6) + W3*V(3,5)
      LEV(1)=1.d0 - (WGHT(1)**2)*SIX1MP*W1*V(1,1)
      LEV(2)=1.d0 - (WGHT(2)**2)*SIX1MP*(W2*V(2,1) + W3*V(2,2))
      TR=LEV(1) + LEV(2)
      DO 20 I=4,NPM1
         W1=1.d0/V(I-2,4)
         W2= -1.d0/V(I-1,4) - 1.d0/V(I-2,4)
         W3=1.d0/V(I-1,4)
         V(I-1,1)=V(I-2,5)*W1 + V(I-2,6)*W2 + V(I-2,7)*W3
         V(I-1,2)=V(I-2,6)*W1 + V(I-1,5)*W2 + V(I-1,6)*W3
         V(I-1,3)=V(I-2,7)*W1 + V(I-1,6)*W2 + V(I,5)*W3
         LEV(I-1)=1.d0 - (WGHT(I-1)**2)*SIX1MP*(W1*V(I-1,1) 
     .             + W2*V(I-1,2) + W3*V(I-1,3))
         TR= TR + LEV(I-1)
   20 CONTINUE
      W1=1.d0/V(NPM2,4)
      W2= -1.d0/V(NPM1,4) - 1.d0/V(NPM2,4)
      W3=1.d0/V(NPM1,4)
      V(NPOINT,1)=V(NPM1,5)*W3
      V(NPM1,1)=V(NPM2,5)*W1 + V(NPM2,6)*W2
      V(NPM1,2)=V(NPM2,6)*W1 + V(NPM1,5)*W2
      LEV(NPM1)=1.d0 - (WGHT(NPM1)**2)*SIX1MP*(W1*V(NPM1,1)
     .             + W2*V(NPM1,2))
      LEV(NPOINT)=1.d0 - (WGHT(NPOINT)**2)*SIX1MP*W3*V(NPOINT,1)
      TR= TR + LEV(NPM1) + LEV(NPOINT)
      RETURN 
      END
