/*
fields, Tools for spatial data
Copyright 2004-2007, Institute for Mathematics Applied Geosciences
University Corporation for Atmospheric Research
Licensed under the GPL -- www.gpl.org/licenses/gpl.html
*/
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <Rmath.h>
#include <float.h>
SEXP addToDiagC(SEXP mat, SEXP numsToAdd, SEXP n)
{
  int In, i;
  double *cMat, *addVals;
  
  //cast R variables to C variables
  In = INTEGER(n)[0];
  cMat = REAL(mat);
  addVals = REAL(numsToAdd);
  
  //Add number to diagonal
  for(i = 0; i < In; i++) {
    cMat[i*In+i] = cMat[i*In+i] + addVals[i];
  }
  
  return R_NilValue;
}
