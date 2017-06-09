/*
c****  # fields is a package for analysis of spatial data written for
c****  # the R software environment .
c****  # Copyright (C) 2017
c****  # University Corporation for Atmospheric Research (UCAR)
c****  # Contact: Douglas Nychka, nychka@ucar.edu,
c****  # National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
c****  #
c****  # This program is free software; you can redistribute it and/or modify
c****  # it under the terms of the GNU General Public License as published by
c****  # the Free Software Foundation; either version 2 of the License, or
c****  # (at your option) any later version.
c****  # This program is distributed in the hope that it will be useful,
c****  # but WITHOUT ANY WARRANTY; without even the implied warranty of
c****  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c****  # GNU General Public License for more details.
*/
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Arith.h>
#include <Rmath.h>
#include <float.h>
SEXP compactToMatCOLD(SEXP compactMat, SEXP len, SEXP n, SEXP diagVal, SEXP lowerTri, SEXP upperTri)
{
  int In, lTri, uTri, i, j, index;
  double dVal;
  double *cMat, *cans;
  
  //cast R variables to C variables
  In = INTEGER(n)[0];
  lTri = INTEGER(lowerTri)[0];
  uTri = INTEGER(upperTri)[0];
  dVal = REAL(diagVal)[0];
  cMat = REAL(compactMat);
  SEXP ans = PROTECT(allocMatrix(REALSXP, In, In));
  cans = REAL(ans);
  
  //set upper or lower triangle of output matrix
  index = 0;
  if(lTri) {
    for(i = 0; i < In; i++) {
      for(j=i+1; j < In; j++) {
        cans[i*In+j] = cMat[index];
        index++;
      }
    }
  }
  index = 0;
  if(uTri) {
    for(i = 0; i < In; i++) {
      for(j=i+1; j < In; j++) {
        cans[j*In+i] = cMat[index];
        index++;
      }
    }
  }
  
  //set diagonal values of output matrix
  for(i = 0; i < In; i++) {
    cans[i*In + i] = dVal;
  }
  
  UNPROTECT(1);
  return ans;
}
