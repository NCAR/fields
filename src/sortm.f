c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
 
  
   
       SUBROUTINE SORTM(K,ki,N)  
C  HEAPSORT ALGORITHM FOR SORTING ON VECTOR OF KEYS K OF LENGTH N      
C  J F MONAHAN        TRANSCRIBED FROM KNUTH, VOL 2, PP 146-7.        
C integer array ki is permuted along with K  
   
      REAL*8 K(N),KK   
      integer ki(N),kki  
      INTEGER R               
      IF(N.LE.1) RETURN      
      L=N/2+1               
      R=N                  
  2   IF(L.GT.1) GO TO 1  
      KK=K(R)  
      kki= ki(R)                 
      K(R)=K(1)  
      ki(R)=ki(1)    
      R=R-1         
      IF(R.EQ.1) GO TO 9     
      GO TO 3               
  1   L=L-1                
      KK=K(L)  
      kki=ki(L)   
  3   J=L         
  4   I=J         
      J=2*J       
      IF(J-R) 5,6,8     
  5   IF(K(J).LT.K(J+1)) J=J+1     
  6   IF(KK.GT.K(J)) GO TO 8      
  7   K(I)=K(J)  
      ki(I)=ki(J)   
      GO TO 4       
  8   K(I)=KK      
      ki(I)=kki   
      GO TO 2    
  9   K(1)=KK   
      ki(1)=kki     
      RETURN       
      END         
