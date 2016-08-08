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
spatialProcess <- function(x, y,  weights = rep(1, nrow(x)),   Z = NULL,
                  mKrig.args = list( m=2),
                cov.function = "stationary.cov", 
                   	cov.args = list(Covariance = "Matern", smoothness = 1),
              	 theta.start = NULL,
                lambda.start = .5, 
                       na.rm = TRUE,
                  	 verbose = FALSE,
                           ...) {
 
# NOTE all ... information is assumed to be for the cov.args list
# overwrite the default choices (some R arcania!)
  ind<- match( names( cov.args), names(list(...) ) )
  cov.args <- c(cov.args[is.na(ind)], list(...))
  if( verbose){
    cat("extra arguments from ... " , names( list(...)) , fill=TRUE)
    cat(" full list from cov.args: ", names(cov.args) )
  }
# NOTE MLEspatialProcess omits NAs
	MLEInfo <- MLESpatialProcess(x, y, weights = weights, Z=Z, 
	                                mKrig.args = mKrig.args,
	                              cov.function = cov.function, 
	                                  cov.args = cov.args,
	                               theta.start = theta.start, 
	                               theta.range = NULL, 
	                                   	gridN = 30,
                            	 lambda.guess = lambda.start,
                                		verbose = verbose
	                                       	 )

	
	if( verbose){
	  cat("Summary from joint optimization", fill=TRUE)
	  print( MLEInfo$MLEJoint$summary )
	  print( MLEInfo$MLEJoint$pars.MLE)
	}
# now fit spatial model with MLE for theta (range parameter)
# reestimate the other parameters for simplicity to get the complete mKrig object
	obj <- do.call( "mKrig", 
	                c( list(x, y, weights=weights,Z=Z) ,
	                          mKrig.args,
	                          list( na.rm=na.rm),
	             list(           
	              cov.function = cov.function,
	              cov.args = cov.args, 
	              lambda = MLEInfo$MLEJoint$pars.MLE[1],
	              theta  = MLEInfo$MLEJoint$pars.MLE[2]
	             )
	            	)
	)
	obj <- c(obj, list(  MLEInfo = MLEInfo, 
	                   theta.MLE =  MLEInfo$MLEJoint$pars.MLE[2],
	                  lambda.MLE = MLEInfo$MLEJoint$pars.MLE[1])
	        )
# replace call to mKrig with this top level one
  obj$call<- match.call()	
	class(obj) <- c( "spatialProcess","mKrig")
 
	return(obj)
}
