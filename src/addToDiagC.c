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
