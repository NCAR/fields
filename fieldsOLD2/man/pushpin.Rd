%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2018
%# University Corporation for Atmospheric Research (UCAR)
%# Contact: Douglas Nychka, nychka@mines.edu,
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
\name{pushpin}
\alias{pushpin}
\title{ Adds a "push pin" to an existing  3-d plot} 
\description{Adds  to an existing 3-d perspective plot 
a push pin to locate a  specific point.}
\usage{
pushpin( x,y,z,p.out, height=.05,col="black",text=NULL,adj=-.1,cex=1.0,...)
}

\arguments{
  \item{x}{x location}
  \item{y}{y location}
  \item{z}{z location}
  \item{p.out}{Projection information returned by persp}   
  \item{height}{Height of pin in device coordinates (default is about
  5\% of the vertical distance ). }
  \item{col}{Color of pin head.}
  \item{text}{Optional text to go next to pin head.}
  \item{adj}{Position of text relative to pin head.}
  \item{cex}{Character size for pin head and/or text}
  \item{\dots}{Additional graphics arguments that are passed to the text 
                     function.}
}

\details{
See the help(text) for the conventions on
the \code{adj} argument and other options for placing text.  

 }
\author{Doug Nychka}

\seealso{drape.plot,persp}

\examples{
# Dr. R's favorite New  Zealand Volcano!
     data( volcano)
     M<- nrow( volcano)
     N<- ncol( volcano)
     x<- seq( 0,1,,M)
     y<- seq( 0,1,,N)

     drape.plot( x,y,volcano, col=terrain.colors(128))-> pm 

     max( volcano)-> zsummit
     xsummit<- x[ row( volcano)[volcano==zsummit]]
     ysummit<- y[ col( volcano)[volcano==zsummit]]

pushpin( xsummit,ysummit,zsummit,pm, text="Summit")

}

\keyword{hplot}
