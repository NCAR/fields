#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP addToDiagC(SEXP, SEXP, SEXP);
extern SEXP compactToMatC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP ExponentialUpperC(SEXP, SEXP, SEXP);
extern SEXP multebC(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RdistC(SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(css)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ddfind)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(dmaket)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(evlpoly)(void *, void *, void *, void *, void *);
extern void F77_NAME(evlpoly2)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(igpoly)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(inpoly)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(mltdrb)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(multrb)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(multwendlandg)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(radbas)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(rcss)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"addToDiagC",        (DL_FUNC) &addToDiagC,        3},
    {"compactToMatC",     (DL_FUNC) &compactToMatC,     6},
    {"ExponentialUpperC", (DL_FUNC) &ExponentialUpperC, 3},
    {"multebC",           (DL_FUNC) &multebC,           8},
    {"RdistC",            (DL_FUNC) &RdistC,            2},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"css",           (DL_FUNC) &F77_NAME(css),           15},
    {"ddfind",        (DL_FUNC) &F77_NAME(ddfind),        10},
    {"dmaket",        (DL_FUNC) &F77_NAME(dmaket),        12},
    {"evlpoly",       (DL_FUNC) &F77_NAME(evlpoly),        5},
    {"evlpoly2",      (DL_FUNC) &F77_NAME(evlpoly2),       7},
    {"igpoly",        (DL_FUNC) &F77_NAME(igpoly),         8},
    {"inpoly",        (DL_FUNC) &F77_NAME(inpoly),         7},
    {"mltdrb",        (DL_FUNC) &F77_NAME(mltdrb),         9},
    {"multrb",        (DL_FUNC) &F77_NAME(multrb),        10},
    {"multwendlandg", (DL_FUNC) &F77_NAME(multwendlandg),  9},
    {"radbas",        (DL_FUNC) &F77_NAME(radbas),         7},
    {"rcss",          (DL_FUNC) &F77_NAME(rcss),          17},
    {NULL, NULL, 0}
};

void R_init_fields(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
