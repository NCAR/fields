# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2018
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
"plot.spatialProcess" <- function(x, digits = 4, which = 1:4, 
    ...) {
    out <- x
    
###########    plot 1  predicted vs. residua values   
    fitted.values <- predict(out)
    if (any(which == 1)) {
        #temp <- summary(out)
        plot(fitted.values, out$y, ylab = "Y", xlab = " predicted values", 
            bty = "n", ...)
        abline(0, 1, col="red")
        hold <- par("usr")
       title("Observations by predicted values")    

    }
##################### plot 2 residual plot
    std.residuals <- (out$residuals * sqrt(out$weights))/out$tau.MLE
    if (any(which == 2)) {
        plot(fitted.values, std.residuals, ylab = "(STD) residuals", 
            xlab = " predicted values", bty = "n", ...)
        yline(0)
    }
################### plot 3 profile over lambda     
    summary<- out$MLEGridSearch$MLEProfileLambda$summary
    profileLambda<- !is.null( summary)
    if (any(which == 3)& profileLambda ) {
    	mar.old<- par()$mar
    
    	# referring to summary[,2] is fragile -- can be either full or REML
    	par( mar= mar.old + c(0,0,0,2) )
            plot(summary$lambda, summary$lnProfileLike.FULL,  xlab = "lambda", 
                ylab ="log Profile Likelihood(lambda)", type="p", log="x", pch=16,cex=.5,
                 ...)
            
            splineFit <- splint(log(summary$lambda), summary$lnProfileLike.FULL, nx = 500)
            lines(exp(splineFit$x),splineFit$y, lwd = 2, col = "red")
            xline( out$lambda.MLE )
            usr.save <- par()$usr
            usr.save[3:4]<- range( summary[,"GCV" ] )
            par( usr= usr.save, ylog=FALSE)
            points(summary$lambda, summary$GCV,
            lty=2, lwd=2, col="blue")
            axis( side=4, col="blue")
            mtext( side=4, line=2,  "GCV function",cex=.75,
                   col="blue")
            title("Profile likelihood over lambda", 
                cex = 0.6)
            par( mar=mar.old)
    }
    
#################### plot 4 profile over theta (range)
    summary<- out$MLEGridSearch$MLEGrid$summary
    profileTheta<- !is.null( summary)
    if ( any(which == 4) & profileTheta ) {
      summary<- out$MLEGridSearch$MLEGrid$summary
    	plot(summary$theta, summary$lnProfileLike.FULL, pch=16, xlab="theta (range parameter)", ylab="log Profile Likelihood (theta)")
    	title("Profile likelihood for theta \n (range parameter)")
    	xline( out$theta.MLE, lwd=2, col="grey40")
    	xline( out$theta.CI, lwd=4, col="grey70", lty=2)
    	splineFit <- splint(summary$theta, summary$lnProfileLike.FULL, nx = 500)
    	lines(splineFit, lwd = 2, col = "red")
           }
}
