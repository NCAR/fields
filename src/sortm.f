c fields, Tools for spatial data
c Copyright (C) 2017, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
 
  
c OLD routine with line numbers and computed if   
c       SUBROUTINE SORTM0(K,ki,N)  
C  HEAPSORT ALGORITHM FOR SORTING ON VECTOR OF KEYS K OF LENGTH N      
C  J F MONAHAN        TRANSCRIBED FROM KNUTH, VOL 2, PP 146-7.        
C integer array ki is permuted along with K  
   
c      double precision K(N),KK   
c      integer ki(N),kki  
c      INTEGER R               
c      IF(N.LE.1) RETURN      
c      L=N/2+1               
c      R=N                  
c  2   IF(L.GT.1) GO TO 1  
c      KK=K(R)  
c      kki= ki(R)                 
c      K(R)=K(1)  
c      ki(R)=ki(1)    
c      R=R-1         
c      IF(R.EQ.1) GO TO 9     
c      GO TO 3               
c  1   L=L-1                
c      KK=K(L)  
c      kki=ki(L)   
c  3   J=L         
c  4   I=J         
c      J=2*J       
c      IF(J-R) 5,6,8     
c  5   IF(K(J).LT.K(J+1)) J=J+1     
c  6   IF(KK.GT.K(J)) GO TO 8      
c  7   K(I)=K(J)  
c      ki(I)=ki(J)   
c      GO TO 4       
c  8   K(I)=KK      
c      ki(I)=kki   
c      GO TO 2    
c  9   K(1)=KK   
c      ki(1)=kki     
c      RETURN       
c      END
      
       SUBROUTINE SORTM(X,ki,N)  
C  HEAPSORT ALGORITHM FOR SORTING ON VECTOR OF KEYS X OF LENGTH N      
C  J F MONAHAN        TRANSCRIBED FROM KNUTH, VOL 2, PP 146-7.        
C integer array ki is permuted along with X    
      double precision X(N),XX
      integer N
      integer ki(N),kki  
      INTEGER R, L,I, J               
      IF(N.LE.1) RETURN      
      L=N/2+1               
      R=N                  
  2   IF(L.GT.1) GO TO 1  
      XX=X(R)  
      kki= ki(R)                 
      X(R)=X(1)  
      ki(R)=ki(1)    
      R=R-1         
      IF(R.EQ.1) GO TO 9     
      GO TO 3               
  1   L=L-1                
      XX=X(L)  
      kki=ki(L)   
  3   J=L         
  4   I=J         
      J=2*J       
C     IF(J-R) 5,6,8
c <  goes here      
      if( (J-R).LT.0) then
        IF(X(J).LT.X(J+1)) J=J+1
      endif
c <= go to here       
      if(  (J-R).LE.0) then 
        IF(XX.GT.X(J)) GO TO 8      
        X(I)=X(J)  
        ki(I)=ki(J)   
      GO TO 4
      endif
c > goes to here       
  8   X(I)=XX      
      ki(I)=kki   
      GO TO 2    
  9   X(1)=XX   
      ki(1)=kki     
      RETURN       
      END
      
