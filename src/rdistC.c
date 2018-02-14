
/*
c****  # fields is a package for analysis of spatial data written for
c****  # the R software environment .
c****  # Copyright (C) 2018
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
SEXP Rdist1C(SEXP x)
     {
         int nx = nrows(x);
	     int dim = ncols(x);
	     void F77_CALL(rdist1)(int *, double *, int *, double *);
         SEXP ans = PROTECT(allocMatrix(REALSXP, nx, nx));
         double *rx = REAL(x), *rans = REAL(ans);
/*         rdist1_( &dim, rx, &nx, rans);    */
         F77_CALL(rdist1)( &dim, rx, &nx, rans);
         UNPROTECT(1);
return ans; 
}

SEXP RdistC(SEXP x1, SEXP x2)
     {
     int n1 = nrows(x1);
	 int n2 = nrows(x2);
	 int dim = ncols(x1);
	 void F77_CALL(rdist)(int *, double *, int *, double *, int *, double *);
         SEXP ans = PROTECT(allocMatrix(REALSXP, n1, n2));
         double *rx1 = REAL(x1),*rx2 = REAL(x2), *rans = REAL(ans);
    	 F77_CALL(rdist)( &dim, rx1, &n1, rx2, &n2, rans);        
         UNPROTECT(1);
return ans; 
}
