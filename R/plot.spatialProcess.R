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
"plot.spatialProcess" <- function(x, digits = 4, which = 1:4, 
    ...) {
    out <- x
    #
    #   don't do plots 2:4 if a fixed lambda
    #
    fitted.values <- predict(out)
    std.residuals <- (out$residuals * sqrt(out$weights))/out$sigma.MLE
    if (any(which == 1)) {
        temp <- summary(out)
        plot(fitted.values, out$y, ylab = "Y", xlab = " predicted values", 
            bty = "n", ...)
        abline(0, 1)
        hold <- par("usr")
       title("Observations by predicted values")    

    }
    if (any(which == 2)) {
        plot(fitted.values, std.residuals, ylab = "(STD) residuals", 
            xlab = " predicted values", bty = "n", ...)
        yline(0)
    }
    if (any(which == 3)) {
    	mar.old<- par()$mar
    	summary<- out$MLEInfo$MLEProfileLambda$summary
    	# referring to summary[,2] is fragile -- can be either full or REML
    	par( mar= mar.old + c(0,0,0,2) )
            plot(summary[,"EffDf" ], summary[,"GCV" ],  xlab = "Eff. number of parameters", 
                ylab = "GCV function", type="l",
                 ...)
            xline( summary[which.min(summary[,"GCV" ] ),"EffDf"])
            usr.save <- par()$usr
            usr.save[3:4]<- range( summary[,2 ] )
            par( usr= usr.save, ylog=FALSE)
            lines(summary[,"EffDf" ], summary[,2 ],
            lty=2, lwd=2, col="blue")
            xline( summary[which.max(summary[,2 ] ),"EffDf"],
                   col="blue")
            axis( side=4)
            mtext( side=4, line=2, "log Profile Likelihood(lamdba)",cex=.75,
                   col="blue")
            title("Profile over lambda", 
                cex = 0.6)
            box()
            par( mar=mar.old)
    }
    if (any(which == 4)) {
      summary<- out$MLEInfo$MLEGrid$summary
      thetaGrid<- (out$MLEInfo$MLEGrid$par.grid)$theta
    	plot(thetaGrid,summary[,2], pch=16, xlab="theta (range parameter)", ylab="log Profile Likelihood (theta)")
    	title("Profile likelihood for theta \n (range parameter)")
    	xline( out$theta.MLE, lwd=2, col="grey")
    	lines( sreg(thetaGrid,summary[,2])$predicted, lwd=2, col="red")
           }
}
