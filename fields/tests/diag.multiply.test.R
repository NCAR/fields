 # fields is a package for analysis of spatial data written for
  # the R software environment .
  # Copyright (C) 2018
  # University Corporation for Atmospheric Research (UCAR)
  # Contact: Douglas Nychka, nychka@ucar.edu,
  # National Center for Atmospheric Research,
  # PO Box 3000, Boulder, CO 80307-3000
  #
  # This program is free software; you can redistribute it and/or modify
  # it under the terms of the GNU General Public License as published by
  # the Free Software Foundation; either version 2 of the License, or
  # (at your option) any later version.
  # This program is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  # GNU General Public License for more details.

suppressMessages(library(fields))
options( echo=FALSE)
set.seed( 234)
test.for.zero.flag<- 1
n <- 5
m <- 4
mat <- array(rnorm(n*m),c(n,m))
mat2 <- array(rnorm(n*m),c(m,n))
vec <- rnorm(n)
vec2 <- rnorm(n)

test.for.zero( mat2 %*% mat, mat2%d*%mat, tol=1e-8 )

test.for.zero( (diag(vec)%*% mat), (vec%d*%mat), tol=1e-8 )


test.for.zero( diag(vec)%*% vec2, vec%d*%vec2,tol=1e-8)
cat("All done with testing diag multiply", fill=TRUE)
options(echo=TRUE)
