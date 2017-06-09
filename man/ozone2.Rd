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

\name{ozone2}
\alias{ozone2}
\title{
  Daily 8-hour ozone averages for sites in the Midwest
}
\description{
The response is 8-hour average (surface) ozone ( from 9AM-4PM) measured in
parts per billion (PPB) for 153 sites in the midwestern US over the period
June 3,1987 through August 31, 1987, 89 days. This season of high 
ozone corresponds with a large modeling experiment using the EPA Regional
Oxidant Model.
}
\usage{
data(ozone2)
}
\format{
The data list has components:
<s-args>
<s-arg name="y">
a 89X153 matrix of ozone values. Rows are days and columns are the 
sites. 
</s-arg>
</s-arg name="lon.lat"> 
Site locations in longitude and latitude as a 153X2 table
</s-arg>
<s-arg name="chicago.subset">
Logical vector indicating stations that form teh smaller 
Chicagoland subset. (see FIELDS ozone data set)
</s-arg>
</s-args>
<s-section name="Reference">
Nychka, D., Cox, L., Piegorsch, W. (1998) Case Studies in Environmental 
Statistics Lecture Notes in Statistics, Springer
Verlag, New York 
}
\examples{
data( ozone2)

# pairwise correlation among all stations
# ( See cover.design to continue this example)
cor.mat<- cor( ozone2$y, use="pairwise")

#raw data image for day number 16 
good<- !is.na( ozone2$y[16,])
out<- as.image( ozone2$y[16,good], x=ozone2$lon.lat[good,])
image.plot( out)
}
\keyword{datasets}
% docclass is data
% Converted by Sd2Rd version 1.21.
