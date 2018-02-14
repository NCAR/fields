%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2018
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

\name{world}
\alias{world}
\alias{world.color}
\alias{in.land.grid}
\alias{world.land}
\title{Plot of the world}
\description{
Plots quickly, medium resolution outlines of large land masses. This
is a simple wrapper for the map function from the maps package. 
}
\usage{
world(...)
world.land(  ...)
world.color( ... )
in.land.grid(...)

}
\arguments{
\item{\dots}{Same arguments used by the \code{map} function from the
maps package.}
}
\details{
See the longstanding \code{maps} package for documentation on this
function.  The functions world.land, world.color and in.land.grid have
been depreciated but can be recovered from versions of fields 6.7.1 or
older.
 }
\seealso{US, in.poly, in.poly.grid}
\examples{
\dontrun{
world()
# add the US
US( add=TRUE,col="blue")

world( fill=TRUE) # land filled in black 

## Western Europe
world( xlim=c(-10,18),ylim=c(36,60),fill=TRUE, col="darkgreen",
border="green1")
}
}
\keyword{hplot}
% docclass is function
