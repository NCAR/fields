/*
c****  # fields is a package for analysis of spatial data written for
c****  # the R software environment .
c****  # Copyright (C) 2018
c****  # University Corporation for Atmospheric Research (UCAR)
c****  # Contact: Douglas Nychka, nychka@ucar.edu,
c****  # National Center for Atmospheric Research,
c****  #        PO Box 3000, Boulder, CO 80307-3000
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
