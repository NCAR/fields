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
spatialProcess <- function(x, y,  weights = rep(1, nrow(x)), 
                           cov.function = "stationary.cov", 
   	cov.args = list(Covariance = "Matern", smoothness = 1),
    	Z = NULL,
	 theta.start=NULL, lambda.start=.5, na.rm=TRUE, verbose=FALSE, ...) {
  if( verbose){
    cat("extra arguments from ... " , names( list(...)) , fill=TRUE) 
  }
# NOTE MLEspatialProcess omits NAs
	MLEInfo <- MLESpatialProcess(x, y, weights=weights, Z=Z, cov.function = cov.function, 
		cov.args = cov.args, theta.start=theta.start, lambda.start = lambda.start,
		verbose=verbose,
		 ...)
	if( verbose){
	  print( MLEfit$summary )
	  print( MLEfit$pars.MLE)
	}
# now fit spatial model with MLE for theta (range parameter)
# reestimate the other parameters for simplicity to get the complete mKrig object
	obj <- mKrig(x, y, weights=weights,Z=Z ,
	              cov.function = cov.function, cov.args = cov.args, 
	            	theta  = MLEInfo$MLEJoint$pars.MLE[2],
	              lambda = MLEInfo$MLEJoint$pars.MLE[1],
	               na.rm = na.rm,
	            	...)
	obj <- c(obj, list(MLEInfo = MLEInfo))
# replace call with this top level one
  obj$call<- match.call()	
	class(obj) <- c( "spatialProcess","Krig")
 
	return(obj)
}
