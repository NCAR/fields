## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  # fig.width = 6, fig.height = 6,
  fig.align = "center"
)
options(digits = 3)

## ----loadpkgs, echo=FALSE, eval=TRUE, message=FALSE---------------------------
library("spam")

## ----echo=TRUE, eval=FALSE, message=FALSE-------------------------------------
#  install.packages("spam")
#  
#  library("spam")

## ----trivial------------------------------------------------------------------
Fmat <- matrix(c(3, 0, 1, 0, 2, 0, 1, 0, 3), nrow = 3, ncol = 3)
Smat <- as.spam(Fmat)

## ----operations---------------------------------------------------------------
Fmat
Smat
Smat %*% t(Smat)
Fmat %*% t(Smat)

## ----nonspam------------------------------------------------------------------
rep(1, 3) %*% Smat

## ----diffbehaviour------------------------------------------------------------
range(Fmat)
range(Smat)

## ----displ1-------------------------------------------------------------------
Smat

## ----displ2-------------------------------------------------------------------
diag.spam(100)

## ----disp3--------------------------------------------------------------------
str(Smat)

## ----disp4--------------------------------------------------------------------
summary(Smat)

## ----disp5, echo=FALSE, warning=FALSE-----------------------------------------
nz <- 2^12
Smat1 <- spam(0, nz, nz)
Smat1[cbind(sample(nz, nz), sample(nz, nz))] <- rnorm(nz)

tmp <- round(summary(Smat1)$density, 3)

## ----disp6, echo=FALSE, warning=FALSE, fig.cap = '\\label{fig:display_spam}Sparsity structure of sparse matrices.'----
par(mfcol=c(1,2),pty='s',mai=c(.8,.8,.2,.2))
display(Smat)
display(Smat1)

## -----------------------------------------------------------------------------
i <- c(2, 4, 4, 5, 5)
j <- c(1, 1, 2, 1, 3)

A <- spam(0, nrow = 5, ncol = 5)
A[cbind(i, j)] <- rep(0.5, length(i))
A <- t(A) + A + diag.spam(5)
A

U <- chol(A)

## ----tree, echo=FALSE, results = 'markup', out.width='40%', fig.show='hold', warning=FALSE, fig.cap='\\label{fig:tree}On the left side the associated graph to the matrix $\\boldsymbol{A}$ is visualized. The nodes of the graph are labeled according to $\\boldsymbol{A}$ (upright) and $\\boldsymbol{P}^T\\boldsymbol{A}\\boldsymbol{P}$ (italics). On the right side the sparsity structure of $\\boldsymbol{A}$ and $\\boldsymbol{P}^T\\boldsymbol{A}\\boldsymbol{P}$ (top row) and the Cholesky factors $\\boldsymbol{R}$ and $\\boldsymbol{U}$ of $\\boldsymbol{A}$ and $\\boldsymbol{P}^T\\boldsymbol{A}\\boldsymbol{P}$ respectively are given in the bottom row. The dashed lines in $\\boldsymbol{U}$ indicate the supernode partition.'----

knitr::include_graphics(c('figures/tree.png', 'figures/ill.png'))

## ----fillin, echo=FALSE, fig.show='hold', warning=FALSE, fig.cap='\\label{fig:ch2:factor}Sparsity structure of the Cholesky factor with MMD, RCM and no permutation of a precision matrix induced by a second-order neighbor structure of the US counties. The values *nnzR* and *fillin* are the number of non-zero elements in the sparsity structure of the factor and the fill-in, respectively.'----

knitr::include_graphics(c('figures/fig_ch2_factors.png'))

