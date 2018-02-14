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
#include <math.h>
#include <R_ext/BLAS.h>
#include <R_ext/Print.h>
#include <unistd.h>

SEXP multebC(SEXP nd, SEXP x1, SEXP n1, SEXP x2, SEXP n2, SEXP par, SEXP c, SEXP work) {
  // NOTE: in original multeb function in fields, "h" was passed but not used and returned as an answer.
  // h is not an argument to this function, but "ans" is equivalent to h and is allocated in C and returned.
  int Ind, In1, In2, r1, r2, d;
  double sum;
  double *Px1, *Px2, *Pc, *Pwork, *cans;
  SEXP ans;
  void expfnC( SEXP, SEXP, SEXP);
  // cast R variables to C variables
  Ind = INTEGER(nd)[0];
  In1 = INTEGER(n1)[0];
  In2 = INTEGER(n2)[0];
  Px1 = REAL(x1);
  Px2 = REAL(x2);
  Pc = REAL(c);
  Pwork = REAL(work);
  // const int tmp = 1;
  const int cn2 = In2;
  
  // allocate answer vector (corresponds to h in fields' multeb)
  ans = PROTECT(allocVector(REALSXP, In1));
  cans = REAL(ans);
  
  // work aray must be dimensioned to size n2
  // outer most loop over columns of x1 and x2 should reduce paging 
  for(r1 = 0; r1 < In1; r1++) {
    // evaluate all basis functions at  x1[r2,.]
    
    for(r2 = 0; r2 < In2; r2++) {
      // zero out sum accumulator
      sum=0.0;
      
      for(d = 0; d < Ind; d++) {
        sum += pow(fabs(Px1[In1*d+r1] - Px2[In2*d+r2]), 2.0);
      }
      
      Pwork[r2]=sum;
    }
    
    // evaluate squared distances with basis functions.
    expfnC(n2, work, par);
    
    // now the dot product you have all been waiting for!
     sum=0.0;
      for(d = 0; d < cn2; d++) {
        sum +=  Pwork[d]*Pc[d] ;
      }
      //  cans[r1] = ddot_(&cn2, Pwork, &tmp, Pc, &tmp);
      cans[r1] = sum;   
  }
  
  UNPROTECT(1);
  return(ans);
}
