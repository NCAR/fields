# spam 2.6-0

SIGNIFICANT USER-VISIBLE CHANGES

* `print.spam()` is able to print non-zero entries of a spam matrix.
* new fast fortran routines used in new `gmult()` to multiply specific spam submatrices with different factors.

INTERNAL CHANGES
* renaming man files such that pkgdown is linking correctly.


# spam 2.5-1

BUG FIXES

* fixing fortran linking warning.


# spam 2.5-0

SIGNIFICANT USER-VISIBLE CHANGES

* improved examples for `spam_random()` and `eigen.spam()`.
* Fortran code used by `eigen.spam()` for non-symmetric matrices is now stable.

INTERNAL CHANGES

* fixing real comparison warnings in Fortran.
* removing remaining and unused print/count/timing variables in Fortran.
* implement dimension upper bound for `eigen.spam()` due to BLAS/LAPACK routines, which are used in Arnoldi iteration (ARPACK).

BUG FIXES

* adjusting spam.Rmd vignette to pandoc2.8.
* correct class checking for matrices.
* cleaning and consolidation of init.c with spam64 (LTO's).


# spam 2.4-0

INTERNAL CHANGES

* cleaning GCC-10 Fortran warnings.


# spam 2.3-0

NEW FEATURES

* New function `spam_random()` to create a random spam matrix.

SIGNIFICANT USER-VISIBLE CHANGES

* Deprecated functions `spam.options()`, `spam.getOption()` are removed.
* Deprecated function `validspamobject()` is now defunct.
* `todo()` and `spam.history()` are removed.
* `summary.spam()` prints whether it is a 32 or 64-bit spam object.

INTERNAL CHANGES

* Fortran modification to address LTO issues.
* Not exported deprecated function `subset.rows.spam()` is now defunct.

BUG FIXES

* Dataset `UScounties.ndorder` contained no-zeros on the diagonal. Now entire diagonal is zero.


# spam 2.2-2

BUG FIXES

* in testthat/test-constructors.R, which uses `base::sample()` (http://developer.r-project.org/blosxom.cgi/R-devel/2019/02/26#n2019-02-26).


# spam 2.2-1

SIGNIFICANT USER-VISIBLE CHANGES

*  New vignette including illustrations and examples
*  Improved documentation for covariance functions like `cov.exp()`.

BUG FIXES

*  `det(spam(1))` bug fix & default of `diag.spam(x)` in `spam_diag()` removed.

INTERNAL CHANGES

*  Spam fit for pkgdown website.
*  Replacing some internal functions by their primitive equivalent.


# spam 2.2-0

SIGNIFICANT USER-VISIBLE CHANGES

*  Implementation of 'eigen.spam()' and 'eigen_approx()' to calculate eigenvalues and eigenvectors for sparse matrices.


# spam 2.1-4

SIGNIFICANT USER-VISIBLE CHANGES

*  'germany.plot' has new argument 'cex.main'.

BUG FIXES

*  'NAOK=TRUE' did not properly dispatch in all matrix multiplications.
*  '.newSpam()' was not always properly called, causing possible errors when creating empty 'spam' objects.

INTERNAL CHANGES

*  Code cleaning and coherence improvements.


# spam 2.1-2

SIGNIFICANT USER-VISIBLE CHANGES

*  In case 64bits are required, issues an error to load package 'spam64' first.
*  Minor fixes in help files and other documentation.

INTERNAL CHANGES

*  Improved registering of compiled functions.
*  Compatibility with new versions of package 'testthat'.


# spam 2.1-1

SIGNIFICANT USER-VISIBLE CHANGES

*  With the addon of the package spam64, we have truly 64 bits available!
*  Dependency on dotCall64.
*  Different option handling, we use now the classical R setting.
    For the moment, the old functionality throws a message but in upcoming versions we use `.Deprecated`.
    'spam.options' -> 'options', 'spam.getOption' -> 'getOption'
*  'validate_spam' superseeds 'validspamobject'.
*  Added the list of Fortran contributors.

INTERNAL CHANGES

*  Many... many. More still to come.
*  Unit testing with testthat.
*  Change of archaic trig functions to standard.
    Minor edits to great circle dist.


# spam 1.4-0

INTERNAL CHANGES

*  Fortran modification to address more pedantic compilations.


# spam 1.3-0

BUG FIXES

*  Adding nonsparse to sparse matrices did not always dispatch properly and caused some errors.
    Thanks to Johan Lindström for pointing out.

SIGNIFICANT USER-VISIBLE CHANGES

*  Help improvements
*  Added additional method dispatching for 'all.equal'.


# spam 1.2-0 and 1.2-1

SIGNIFICANT USER-VISIBLE CHANGES

*  Renaming of demos from JSS10 paper. 

*  Continuing to implement 'spam_xxx' function names: Like 'spam_rdist' there is now 'spam_diag'.
*  'subset.rows.spam' deprecated.
    There is an internal function 'subset_rows.spam'.

INTERNAL CHANGES

*  Ample minor modifications in 'DESCRIPTION' and 'NAMESPACE' for CRAN check conformity.
*  Elimination of a few unnecessary functions (not exported).
*  Minor modifications in tests.


# spam 1.1-0

SIGNIFICANT USER-VISIBLE CHANGES

*  Upload of vignette linked to JSS15 paper. Related demos and tests.
*  "DUP=FALSE" is deprecated and will be disabled in future versions of R.
    All "DUP" arguments have been eliminated.
    You have to expect a slightly slower version of spam.


# spam 1.0 and 1.0-1

*  This version is up to 'DESCRIPTION' and its implied changes identical to 0.90-1. With the upcoming JSS article "Pitfalls in the implementation of Bayesian hierarchical modeling of areal count data.
    An illustration using BYM and Leroux models.
    " a "major" version jump is adequate.
*  Referencing to spam data through spam::...
*  Implemented `1.1.3.1 Suggested packages` approach.


# spam 0.70/0.80/0.90/0.90-1

SIGNIFICANT USER-VISIBLE CHANGES

*  Introduction of many as('spam','...') functions.
*  Coercion function `as.vector` for spam objects.
*  Wrapper functions `spam_rdist` and 'spam_rdist.earth` for smooth use in `fields`.
*  The use of `update(A, B)` without assignment has been eliminated.
    This is one way to address the change in memory handling changes from R 3.0.2 to R 3.1.0.
    There is a slight overhead in memory. If this causes problems, let me know.
*  Adjustment of the license.

NEW FEATURES

*  Arguments `diag` and `eps` in `nearest.dist` cause now an error.
*  Further augmented help pages.

BUG FIXES

*  The demo now points to the new JSS article.

INTERNAL CHANGES

*  Set 'structurebased=TRUE' for the demos.
*  Link to upcoming JSS article in one of the demos.
*  'update.spam.chol.NgPeyton' preserves the structure (pointed out by Chris Paciorek), see above.
*  Using similar License approach as SparseM.
    New files `README`, `inst/0LICENSE`.
*  File renaming (OChangeLog -> 0ChangeLog) 
*  Adjusted error messages for precmat.RW2


# spam 0.60-0

SIGNIFICANT USER-VISIBLE CHANGES

*  Using the flag 'structurebased', the behavior of spam is now more consistent.
*  "Arith", "Compare", "Logic" (getGroupMembers("Ops")) have now a consistent behavior.

NEW FEATURES

*  Few new S3 functions for simplicity: 'var.spam', 'eigen.spam', ...
*  New constructor functions 'colindices<-' etc.
    Maybe additional tests may be required.
*  Operators from 'Arith' obey now the structure based calculation.

BUG FIXES

*  'inefficiencywarning' passes message correctly.

INTERNAL CHANGES

*  many more spam/tests/*.
*  Consistent use of 'spam' and 'vector' siglist for 'Ops'.
*  Minor cleaning of Fortran code.
*  Renaming/restructuring/cleaning of files...
*  Fortran arguments are copied when updating the cholesky structure.


# spam 0.50-0

SIGNIFICANT USER-VISIBLE CHANGES

*  Using the flag 'structurebased=FALSE', the behavior of spam is now much, much closer to regular matrix calculations.
    This is illustrated when calculating gamma of a sparse matrix.
*  Along the same lines, the flag 'NAOK=TRUE' allows the use of the "not finite numbers" (NA, NaN, Inf).
    We have tested many, many functions but full fledged use is not yet guaranteed.
*  Currently, we can still guarantee backwards compatibility...

NEW FEATURES

*  New functions 'crossprod' and 'tcrossprod' as well as according method definitions.
*  New constructor functions 'rowpointers<-' etc.
*  Better option handling. The option 'safemode' is now 'safemodevalidity'. Additionally, new option 'NAOK'.
*  Help pages have been improved.
*  Operators from 'Summary' and 'Math' obey now the structure based calculation.
    ('Math2' inherently does).

BUG FIXES

*  rmvnorm.[].const now work properly for any number of constraints and n.
*  Assignment handles properly recycling.
*  todo() now works properly.

INTERNAL CHANGES

*  eliminated {d,i}check by equivalent coercion.
*  Consistent use of NAOK in Fortran calls.
*  Minor cleaning of Fortran code.
*  Renaming/restructuring of files...


# spam 0.42-0

NEW FEATURES

*  More consistent handling of subsetting. 
    Warning is issued if subsetting with NA

BUG FIXES

*  Fixed several issues when rowsubsetting...

INTERNAL CHANGES

*  Additional tests for positive definiteness in 'chol'.


# spam 0.41-0

NEW FEATURES

*  Functions grid_trace2() received more functionality.

BUG FIXES

*  Eliminated bug in cov.mat().
    Pointed out by Joshua French.

INTERNAL CHANGES

*  Updated DESCRIPTION: added Florian Gerber [ctb].
*  Minor code and help cleanup. Additional testing files.
    File header edits.
*  Addressed Rdevel CMD check --as-cran Notes, especially workaround for DUP=FALSE.


# spam 0.40-0

BUG FIXES

*  A severe bug in subsetting a spam object with a nx2 matrix crept in spam in version 0.29-3.
    Thanks to Andrew Hong and Beat Briner for pointing out.
*  To simplify communication, we have switched increased the tenth version number.
*  All other changes are of cosmetic nature.


# spam 0.30-x

SIGNIFICANT USER-VISIBLE CHANGES

*  Added several plots to visualize several MCMC chains ('grid_trace2', 'grid_zoom', ...).

NEW FEATURES

*  New function 'germany.plot' to draw the landkreise.
    ('map.landkreis' is now obsolete).
*  Switched from 'tim.colors' to 'colorRampPalette' in 'germany.plot'.
*  Metadata in 'germany.info', polygon definitions in 'germany.poly' ('germany' kept for backwards compatibility).

INTERNAL CHANGES

*  Switched to mercurial for maintaining the package.
*  Updated ChangeLog file (hg log).
*  Increased dependency to >= R 2.15.
*  Minor code and help cleanup.


# spam 0.29-0, 0.29-1, 0.29-2, 0.29-3

SIGNIFICANT USER-VISIBLE CHANGES

*  There is a generic conflict with 'backsolve' between spam and other packages (e.g., bdsmatrix).
   To avoid the issue, we use the standard generic implemented in 'methods' which requires an additional argument for version 0.29-0 (see also PR#14883).
    However to maintain backwards compatibility with packages that depend on spam, this was reverted in 0.29-1.
    Currently, this conflict is not properly solved.
    I propose  to load 'spam' first then the other packages, followed by manually calling:

  *    setMethod("backsolve","spam.chol.NgPeyton",backsolve.spam)
        setMethod("backsolve","spam",backsolve.spam)

    Stay tuned...

*  Calls like:
 
  *    mat <- diag.spam(4)
        diag(mat[-1, ]) <- 3
        diag.spam(mat[ , -1]) <-2

    now work. They are, however, somewhat inefficient.
    'toeplitz.spam' is to be prefered.
    Pointed out by Florian Gerber.

*  The Gibbs sampler in the demo article-jss-example2 contains several bugs, pointed out by Steve Geinitz and Andrea Riebler.
    I'll post an updated sampler in a future release.

NEW FEATURES

*  New functions 'rmvnorm.const', 'rmvnorm.prec.const' and 'rmvnorm.canonical.const' to draw constrained multivariate normal variates.
*  New functions 'precmat' (wrapper to), 'precmat.RW1', 'precmat.RW2', 'precmat.season', 'precmat.IGMRFreglat' and 'precmat.IGMRFirreglat' to create precision matrices for IGMRF.
*  New methods 'rowSums', 'colSums' and  'rowMeans', 'colMeans' for 'spam' objects.
*  New methods 'head' and 'tail' for 'spam' and 'spam.chol.NgPeyton' objects.
*  New method 'chol2inv' for 'spam' object.
*  New option 'inefficiencywarning': handling of warnings issued in case of an inefficient calculation. 
*  New option 'structurebased': should operations be performed on the nonzero entries or on including the zeros.
    Classical example: what should the cosine of a sparse matrix look like?
    In the near future, all operations from Math and Ops will include this option.
    Some loss of backwards compatibility might be lost in the future.

INTERNAL CHANGES

*  New much faster approach to extract rows. For not too sparse large matrices improvements over two orders of magnitudes are achieved.
*  Elininated '.Internal()' calls that induce a 'Note' on CRAN checks. This also implied a minor rewrite of 'image.spam'.
*  Minor code improvements.
*  Eliminated non-API calls (29.1).
*  Rewritten .C("bincode",...) call as suggested by Brian Ripley (29.2).

BUG FIXES

*  Bug fix that occures when multiplying elementwise matrices that have non-intersecting structures (pointed out by Corentin Barbu).
*  Bug fix in triangular backsolves involving 'spam' objects and rhs matrices.
*  Bug fix in triangular backsolve causing errors on some architectures.


# spam 0.28

NEW FEATURES

*  New function 'cleanup' (suggested by Simon Barthelme).
*  Extending help files.
*  Improved functionality of 'isSymmetric'.

INTERNAL CHANGES

*  Proper storage of data files.
*  Cleaning up argument names within spam functions.
*  Cleaning up old Fortran code, i.e., eliminating unnecessary subroutines and write calls.

BUG FIXES

*  Bug fix that may occure when extracting zero elements (pointed out by Corentin Barbu).


# spam 0.27

NEW FEATURES

*  Requires now R2.10 and higher.
*  Functions to create Toeplitz and circulant matrices.
*  Function to create precision matrices for gridded GMRF.
*  Improvements in the mle.* functions.
*  Method diff for sparse matrices (suggested by Paul Eilers).
*  Improvement of help pages.
*  Eliminated some help aliases to base functions (for which no 'usage' is given).


INTERNAL CHANGES

*  Change to iL coding.
*  Start to using 'identical'.
*  Code cleaning due to requirement of R2.10 and higher.

BUG FIXES

*  Bug fix in as.spam.list (thanks to Paul Eilers).
*  Bug fix in demo(spam) (thanks to Thomas Gsponer).


# spam 0.24, 0.25 and 0.26

*  Devel versions, not released.


# spam 0.23

NEW FEATURES

*  Further improved versions of demos.
*  Some improvements to meet Rd standards.
    Adjustments for future R versions.


# spam 0.22

NEW FEATURES

*  Improved versions of demos. Synchronized with the JSS article.
*  Additional changes and improvements in the help files (thanks to Steve Geinitz).


# spam 0.21

NEW FEATURES

*  New NEWS file, to work better with news() command.
    The previous is available under ONEWS.
*  New functions bandwidth, permutation, mle[.nomean][.spam], neg2loglikelihood[.spam].
*  Renamed adiag.spam to bdiag.spam.
*  Cleaned up argument naming with the rmvnorm.* suite.

INTERNAL CHANGES

*  Various Fortran code, R code and help file improvements.

BUG FIXES

*  Minor change in one of the demos (solves a 64bit issue).
