/*
fieldsMAGMA
Copyright 2004-2015, Institute for Mathematics Applied Geosciences
University Corporation for Atmospheric Research
Licensed under the GPL -- www.gpl.org/licenses/gpl.html
*/
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <Rmath.h>
#include <float.h>

SEXP ExponentialUpperC(SEXP distMat, SEXP n, SEXP alpha)
{
  int In, i, j;
  double dAlpha;
  double *dMat, *cans;
  
  //cast R variables to C variables
  In = INTEGER(n)[0];
  dAlpha = REAL(alpha)[0];
  dMat = REAL(distMat);
  SEXP ans = PROTECT(allocMatrix(REALSXP, In, In));
  cans = REAL(ans);
  
  //intialize entire array to zero DWN May-4-2016
  for(i = 0; i < (In*In); i++) {
        cans[i]= 0.0;
    }
  
  //set upper triangle of output matrix
  for(i = 0; i < In; i++) {
    for(j=0; j<= i; j++) {
      if(i == j)
        cans[i*In+j] = 1.0;
      else
        cans[i*In+j] = exp(-1*dMat[i*In+j]*dAlpha);
    }
  }
  
  UNPROTECT(1);
  return ans;
}
