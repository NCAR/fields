# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    

mKrig.MLE.joint <- function(x, y, weights = rep(1, nrow(x)), 
                            lambda.guess = 1, cov.params.guess=NULL, 
                            cov.fun="stationary.cov", cov.args=NULL, 
                            Z = NULL, optim.args=NULL, find.trA.MLE = FALSE, 
                            ..., verbose = FALSE) {
  
  warning(paste0("mKrig.MLE.joint is deprecated and may be removed in a ", 
                 "future release.  Use mKrigMLEJoint instead."))
  
  do.call("mKrigMLEJoint", c(list(x,y, weights, lambda.guess, cov.params.guess, 
                                  cov.fun, cov.args, Z, optim.args, find.trA.MLE, 
                                  verbose=verbose), list(...)))
}