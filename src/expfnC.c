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

SEXP expfnC(SEXP n, SEXP d2, SEXP par) {
  int In, i;
  double Dpar, par2;
  double *Pd2;
  
  //caste R variables to C variables, allocate answer vector
  In = INTEGER(n)[0];
  Dpar = REAL(par)[0];
  par2 = Dpar/2;
  Pd2 = REAL(d2);
  
  for(i = 0; i < In; i++) {
    Pd2[i] = exp(-1*pow(Pd2[i], par2));
  }
  
  return(R_NilValue);
}
