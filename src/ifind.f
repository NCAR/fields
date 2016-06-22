c fields, Tools for spatial data
c Copyright 2015, Institute for Mathematics Applied Geosciences
c University Corporation for Atmospheric Research
c Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
  
      INTEGER FUNCTION IFIND(X,XK,N)                                   
C  FIND I SUCH THAT XK(I) LE X LT XK(I+1)                              
C  IFIND=0  IF X LT XK(1)                                              
C  IFIND=N  IF X GT XK(N)                                              
C  J F MONAHAN  JAN 1982  DEPT OF STAT, N C S U, RALEIGH, NC 27650     
      REAL*8 X,XK(N)                                                   
      IFIND=0                                                          
      IF(X.LT.XK(1)) RETURN                                            
      IFIND=N                                                          
      IF(X.GE.XK(N)) RETURN                                            
      IL=1                                                             
      IU=N                                                             
  1   IF(IU-IL.LE.1) GO TO 4                                           
      I=(IU+IL)/2                                                      
C      IF(X-XK(I)) 2,5,3                                                
      IF( (X-XK(I)).eq.0) go to 5
      IF( (X-XK(I)).gt.0)  go to 3                                                
  2   IU=I                                                             
      GO TO 1                                                          
  3   IL=I                                                             
      GO TO 1                                                          
  4   IFIND=IL                                                         
      RETURN                                                           
  5   IFIND=I                                                          
      RETURN                                                           
      END                                                              
