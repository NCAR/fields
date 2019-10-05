
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <Rmath.h>
#include <float.h>
#include <math.h>


#define min(a,b) (a<b?a:b)

SEXP distMatHaversin(SEXP p1, SEXP radius, SEXP ans)
{
  // Florian Gerber, gerber@mines.edu
  // Calculates the great-cricle distance matrix for one set of locations.
  // The formula from the R function geosphere::distHaversin() v1.5-10 is used.
  // Note:
  // The memory for 'ans' has to be allocated in R and the vector
  // should be initialized with zeros.
  // If this is done in C with
  //      PROTECT(ans = allocMatrix(REALSXP, n1, n1));
  //      ...
  //      UNPROTECT(1); 
  //      return(ans);
  // the memory usage is doubled (why?). 
  
  R_len_t i, j, n1 = length(p1) / 2;
  double *rp1 = REAL(p1);
  double *rradius=REAL(radius), *rans=REAL(ans); // 'ans': vector of zeros
  double p1Lon, p1Lat, p2Lon, p2Lat, dLon, dLat, tmp;
  const double pi = 3.141592653589793238462643383280;
  const double torad = pi / 180.0;  

  for(i = 0; i < n1; i++) {
    for(j = i+1; j < n1; j++){
      p1Lon = rp1[i] * torad;
      p1Lat = rp1[i+n1] * torad;
      p2Lon = rp1[j] * torad;
      p2Lat = rp1[j+n1] * torad;
      dLon = (p1Lon - p2Lon) / 2.0;
      dLat = (p1Lat - p2Lat) / 2.0;
      tmp = pow(sin(dLat), 2) + cos(p1Lat) * cos(p2Lat) * pow(sin(dLon), 2);
      tmp = min(tmp, 1.0);
      tmp = 2.0 * atan2(sqrt(tmp), sqrt(1.0 - tmp)) * rradius[0];
      rans[i + n1*j] = tmp;
      rans[j + n1*i] = tmp; //matrix is symmetric
    }
  }
  return(R_NilValue);
}


SEXP distMatHaversin2(SEXP p1, SEXP p2, SEXP R, SEXP ans)
{
  // Florian Gerber, gerber@mines.edu
  // Calculates the great-cricle distance matrix between two sets of locations.
  // The formula from the R function geosphere::distHaversin() v1.5-10 is used.

  R_len_t i, j, n1 = length(p1) / 2, n2 = length(p2) / 2;
  double *rp1 = REAL(p1), *rp2 = REAL(p2);
  double *rR=REAL(R), *rans=REAL(ans);
  double p1Lon, p1Lat, p2Lon, p2Lat, dLon, dLat, tmp;
  const double pi = 3.141592653589793238462643383280;
  const double torad = pi / 180.0;
  
  for(i = 0; i < n1; i++) {
    for(j = 0; j < n2; j++){
      p1Lon = rp1[i] * torad;
      p1Lat = rp1[i+n1] * torad;
      p2Lon = rp2[j] * torad;
      p2Lat = rp2[j+n2] * torad;
      dLon = (p1Lon - p2Lon) / 2.0;
      dLat = (p1Lat - p2Lat) / 2.0;
      tmp = pow(sin(dLat), 2) + cos(p1Lat) * cos(p2Lat) * pow(sin(dLon), 2);
      tmp = min(tmp, 1.0);
      tmp = 2.0 * atan2(sqrt(tmp), sqrt(1.0 - tmp)) * rR[0];
      rans[i + n1*j] = tmp;
    }
  }
  return(R_NilValue);
}
