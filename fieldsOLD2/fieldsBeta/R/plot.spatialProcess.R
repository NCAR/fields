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
    tauMLE<- out$summary["tau"]
    std.residuals <- (out$residuals * sqrt(out$weights))/tauMLE
    if (any(which == 2)) {
        plot(fitted.values, std.residuals, ylab = "(STD) residuals", 
            xlab = " predicted values", bty = "n", ...)
        yline(0)
    }
################### plot 3 profile over lambda 
    profileLambda<- !is.null( out$lambdaProfile)
    
    if (any(which == 3)& profileLambda ) {
      summary<- out$lambdaProfile$summary
    	mar.old<- par()$mar
    
    	# referring to summary[,2] is fragile -- can be either full or REML
    	par( mar= mar.old + c(0,0,0,2) )
            plot(summary$lambda, summary$lnProfileLike.FULL,  xlab = "lambda", 
                ylab ="log Profile Likelihood(lambda)",
                type = "p", log="x", pch=16,cex=.5,
                 ...)
            
            splineFit <- splint(log(summary$lambda), summary$lnProfileLike.FULL, nx = 500)
            lines(exp(splineFit$x),splineFit$y, lwd = 2, col = "red")
            xline( out$lambda.MLE )
            usr.save <- par()$usr
            RGCV<-range( summary[,"GCV" ] )
# 8% expansion of scale 
            delta<- (RGCV[2]- RGCV[1])* .08
            usr.save[3:4]<-  c( RGCV[1] -delta, RGCV[2]+ delta)
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
    
#################### plot 4 profile over aRange (range)
    profileARange<- !is.null(out$aRangeProfile)
    if ( any(which == 4) & profileARange ) {
      summary<- out$aRangeProfile$summary
    	plot(summary$aRange, summary$lnProfileLike.FULL, pch=16, xlab="aRange (range parameter)", ylab="log Profile Likelihood (aRange)")
    	title("Profile likelihood for aRange \n (range parameter)")
    	xline( out$aRange.MLE, lwd=2, col="grey40")
    #	xline( out$aRange.CI, lwd=4, col="grey70", lty=2)
    	splineFit <- splint(summary$aRange, summary$lnProfileLike.FULL, nx = 500)
    	lines(splineFit, lwd = 2, col = "red")
           }
}
