%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2016
%# University Corporation for Atmospheric Research (UCAR)
%# Contact: Douglas Nychka, nychka@ucar.edu,
%# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
%#
%# This program is free software; you can redistribute it and/or modify
%# it under the terms of the GNU General Public License as published by
%# the Free Software Foundation; either version 2 of the License, or
%# (at your option) any later version.
%# This program is distributed in the hope that it will be useful,
%# but WITHOUT ANY WARRANTY; without even the implied warranty of
%# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%# GNU General Public License for more details.
%#
%# You should have received a copy of the GNU General Public License
%# along with the R software environment if not, write to the Free Software
%# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%# or see http://www.r-project.org/Licenses/GPL-2    

\name{arrow.plot}
\alias{arrow.plot}
\title{
  Adds arrows to a plot
}
\description{
Adds arrows at specified points where the arrow lengths are scaled to
fit on 
the plot in a reasonable manner. A classic use of this function is to
depict a vector field. At each point (x,y) we have a vector with
components (u,v).  Like the arrows function this adds arrows to an 
existing plot.  
}
\usage{
arrow.plot(a1, a2, u = NA, v = NA, arrow.ex = 0.05,  
     xpd = TRUE, true.angle = FALSE, arrowfun=arrows,...)
}
\arguments{
\item{a1}{
The x locations of the tails of the arrows or a  2 column matrix giving 
the x and y coordinates of the arrow tails. 
}
\item{a2}{
The y locations of the tails of the arrows or a 2 column matrix giving  
the u and v coordinates of the arrows. 
}
\item{u}{
The u components of the direction vectors if they are not specified in 
the a1 argument 
}
\item{v}{
The v  components of the direction vectors if they are not specified in     
the a2 argument 
}
\item{arrow.ex}{
Controls the length of the arrows. The length is in terms of the fraction 
of the shorter axis in the plot. So with a default of .05 20 arrows  
of maximum length can line up end to end along the shorter axis. 
}
\item{xpd}{
If true does not clip arrows to fit inside the plot region, default is 
not to clip. 
}
\item{true.angle}{
If true preserves the true angle of the (u,v) pair on the plot. E.g. if  
(u,v)=(1,1) then the arrow will be drawn at 45 degrees.  
}
\item{arrowfun}{
The actual arrow function to use. The default is standard R \code{arrows}.
However,  Tamas K Papp suggests  \code{p.arrows} from sfsmisc which makes prettier 
arrows.
}


\item{\dots}{
Graphics arguments passed to the arrows function that can can change the 
color or arrow sizes. See help on this for 
details. 
}
}
\details{
This function is useful because (u,v) may be in very different scales
from the locations (x,y). So some careful scaling is needed to plot the
arrows. 
The only tricky thing about this function is whether you want the true 
angles on the plot. For overlaying a  vector field on top of contours 
that are the streamlines  true.angle should be false. In this case you 
want u and v to be scaled in the same way as the x and y variables.  
If the scaling is not the same then the arrows will not look like tangent 
vectors to the streamlines.  
An application where the absolute angles are meaningful might be the hands of a 
clock showing different times zones on a world map. Here true.angle=T is
appropriate, the clock hands should preserve the right angles. 
}
\seealso{arrows} 

\examples{
#
# 20 random directions at 20 random points

x<- runif( 20)
y<- runif( 20)
u<- rnorm( 20)
v<- rnorm( 20)
plot( x,y)
arrow.plot( x,y,u,v) # a default that is unattractive 

plot( x,y, type="n")
arrow.plot( x,y,u,v, arrow.ex=.2, length=.1, col='green', lwd=2) 
# thicker lines in green, smaller heads and longer tails. Note length, col and lwd are
# options that the arrows function itself knows about. 
  
}
\keyword{aplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
