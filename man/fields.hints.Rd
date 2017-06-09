%# fields  is a package for analysis of spatial data written for
%# the R software environment .
%# Copyright (C) 2017
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

\name{fields.hints}
\alias{fields.hints}
\alias{fields.style}
\alias{fields.color.picker}
\title{
fields - graphics hints  
}

\description{
 Here are some technical hints for assembling multiple plots with common
legends or axes and setting the graphics parameters to make more
readable figures. Also we an index to the defaultcolors in R graphics
and setting their definitions in LaTeX.
All these hints use the standard graphics environment. 
}
\usage{
fields.style()
fields.color.picker()
}
\details{
\code{fields.style} is a simple 
function to enlarge the characters in a plot and set the colors. List this out to modify the choices.

\preformatted{

##Two examples of concentrating a panel of plots together
## to conserve the white space. 
## see also the example in image.plot using split.screen. 
## The basic trick is to use the oma option to reserve some space around the 
## plots.  Then unset the outer margins to use that room. 

library( fields)

# some hokey image data
x<- 1:20
y<- 1:15
z<-  outer( x,y,"+")
zr<- range( c(z))

# add common legend to 3X2 panel 

par( oma=c(4,0,0,0))
set.panel( 3,2)
par( mar=c(1,1,0,0))

# squish plots together with just 1 space between
for ( k in 1:6){
image( x,y,z, axes=FALSE, xlab="", ylab="", zlim=zr)
}

par( oma=c(0,0,0,0))
image.plot( zlim=zr, legend.only=TRUE, horizontal=TRUE, legend.mar=5)

# you may have to play around with legend.mar and the oma settings to
# get enough space.



##
### also add some axes on the sides. in a lattice style 
## note oma adds some more room at bottom. 

par( oma=c(8,6,1,1))
set.panel( 3,2)
par( mar=c(1,1,0,0))
##
for ( k in 1:6){
 image( x,y,z, axes=FALSE, xlab="", ylab="", zlim=zr)
 box() # box around figure

# maybe draw an x axis 
  if( k  \%in\% c(5,6) ){
  axis( 1, cex.axis=1.5)
  mtext( line=4, side=1, "Xstuff")}

# maybe draw a y axis
  if( k  \%in\% c(1,3,5) ){
  axis( 2, cex.axis=1.5)
  mtext( line=4, side=2, "Ystuff")}
}

# same trick of adding a legend strip. 
par( oma=c(0,0,0,0))
image.plot( zlim=zr, legend.only=TRUE, horizontal=TRUE, legend.mar=5)

# reset panel 
set.panel()


####
# show colors
## the factory colors:

clab<- colors()
n<- length( clab)
N<- ceiling( sqrt(n) )
M<- N
temp<- rep( NA,M*N)
temp[1:n] <- 1:n
z<- matrix(temp, M,N)

image(seq(.5,M+.5,,M+1), seq(.5,N+.5,,N+1)
       , z,  col=clab, axes=FALSE, xlab="", ylab="")


# see the  function  fields.color.picker() to locate colors



# dumping out colors by name for a latex document
# this creates text strings that are the LaTeX color definitions
# using the definecolor function. 

# grab all of the R default colors 
clab<- colors()

for( nn in clab){
  temp<- signif(col2rgb(nn)/256, 3)
   cat(
    "\\definecolor{", 
                nn, "}",
    "{rgb}{", temp[1], 
          ",", temp[2], 
          ",", temp[3], 
           "}", fill=TRUE , sep="")
 }

# this loop prints out definitions such as 
# \definecolor{yellowgreen}{rgb}{0.602,0.801,0.195}
# having loaded the color package in LaTeX 
# defining this color 
# use the construction  {\color{yellowgreen}  THIS IS A COLOR}
# to use this color in a talk or document. 

# this loop prints out all the colors in LaTeX language
# as their names and can be converted to a pdf for handy reference.

sink( "showcolors.tex")

clab<- colors()
for( nn in clab){
  temp<- signif(col2rgb(nn)/256, 3)
   cat(
    "\\definecolor{", 
                nn, "}",
    "{rgb}{", temp[1], 
          ",", temp[2], 
          ",", temp[3], 
           "}", fill=TRUE , sep="")
   cat( paste("{ \\color{",nn,"} ", nn," $\\bullet$ \\\\ }", sep=""),
                    fill=TRUE)
}
sink()

} %end preformatted

} 

\keyword{hplot}
