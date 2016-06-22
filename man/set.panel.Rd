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

\name{set.panel}
\alias{set.panel}
\title{
  Specify a panel of plots 
}
\description{
Divides up the graphics window into a matrix of plots. 
}
\usage{
set.panel(m=1, n=1, relax=FALSE)
}
\arguments{
\item{m}{
Number of rows in the panel of plots 
}
\item{n}{
Number of columns in the panel. 
}
\item{relax}{
If true and the par command is already set for multiple plots,  
then the set.panel command is ignored. The default is  
relax set to false.  
}
}
\details{
After set.panel is called, the graphics screen is reset to 
put  plots  according to a m x n table. Plotting starts in 
the upper left hand corner and proceeds row by row.  After 
m x n plots have been drawn, the next plot will erase the 
window and start in the 1,1 position again.  This function 
is just a repackaging for specifying the mfrow argument to 
par. 
Setting up a panel of plots is a quick way to  change  the 
aspect  ratio  of  the graph (ratio of height to width) or 
the size. For example, plotting 2 plots to a page produces 
a  useful  size  graph  for including in a report. You can 
print out the graphs at any stage without having  to  fill 
up the entire window with plots. This function, except for the "relax"
option is equivalent to the  S sequence: par( mfrow=c(m,n)).
}
\section{Side Effects}{
The function will echo your choice of m and n to the terminal.  
}
\seealso{
\code{par}  
}
\examples{
set.panel(5,2) #divide screen to hold 10 plots where there are 5 rows  
	       #and 2 columns 
plot( 1:10) 
plot( 2:8)

set.panel() #reset screen to one plot per screen 
}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
