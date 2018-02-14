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
\name{WorldBankCO2}
\alias{WorldBankCO2}
\docType{data}
\title{Carbon emissions and demographic covariables by country for 1999.}
\description{
These data are a small subset of the demographic data compiled by the World Bank. The
data has been restricted to 1999 and to countries with a population larger than 1 million. Also, only countries reporting all the covariables are included. 
}
 \usage{ data(WorldBankCO2)}
\format{
  This a 75X5 matrix with the row names identifying  countries  and 
columns the covariables: 
\code{ "GDP.cap" "Pop.mid" "Pop.urb" "CO2.cap" "Pop"}
\itemize{
\item GDP.cap: Gross domestic product (in US dollars) per capita.
\item Pop.mid: percentage of the population within the ages of 15 through 65.
\item Pop.urb: Precentage of the population living in an urban environment
\item CO2.cap: Equivalent CO2 emmissions per capita 
\item Pop: Population
}

}
\section{Reference}{ 
Romero-Lankao, P., J. L. Tribbia and D. Nychka (2008) Development and greenhouse 
gas emissions deviate from the modernization theory and convergence hypothesis. Cli- 
mate Research 38, 17-29. 
}
\examples{
data(WorldBankCO2)
plot( WorldBankCO2[,"GDP.cap"], WorldBankCO2[,"CO2.cap"], log="xy")
}
\section{Creating dataset}{
Listed below are scripts to create this data set from spread sheet on the
World Bank CDs:
\preformatted{
## read in comma delimited spread sheet
  read.csv("climatedemo.csv", stringsAsFactors=FALSE)->hold
## convert numbers to matrix of data
  Ddata<- as.matrix(  hold[,5:51] )
  Ddata[Ddata==".."] <- NA
## still in character form parse as numeric
  Ddata<- matrix( as.numeric( Ddata), nrow=1248, ncol=ncol( Ddata),
  dimnames=list( NULL, format( 1960:2006) ))
## these are the factors indicating the different variables
### unique( Fac) gives the names of factors 
  Fac<- as.character( hold[,1])
  years<- 1960:2006
# create separate tables of data for each factor
 temp<- unique( Fac)
## also subset Country id and name
  Country.id<- as.character( hold[Fac== temp[1],3])
  Country<- as.character( hold[Fac== temp[1],4])
  Pop<- Ddata[ Fac== temp[2],]
  CO2<- Ddata[ Fac== temp[1],]
  Pop.mid<- Ddata[ Fac== temp[3],]
  GDP.cap<- Ddata[ Fac== temp[4],]
  Pop.urb<- Ddata[ Fac== temp[5],]
  CO2.cap<- CO2/Pop
  dimnames( Pop)<- list( Country.id,format(years))
  dimnames( CO2)<- list( Country.id,format(years))
  dimnames( Pop.mid)<- list( Country.id,format(years))
  dimnames( Pop.urb)<- list( Country.id,format(years))
  dimnames( CO2.cap)<- list( Country.id,format(years))
# delete temp data sets
  rm( temp)
  rm( hold)
  rm( Fac)
# define year to do clustering.
  yr<- "1999"
# variables for clustering combined as columns in a matrix
  temp<-cbind( GDP.cap[,yr], Pop.mid[,yr], Pop.urb[,yr],CO2[,yr],Pop[,yr])
# add column names and figure how many good data rows there are.
        dimnames( temp)<-list( Country, c("GDP.cap","Pop.mid","Pop.urb",
                                          "CO2.cap", "Pop"))
        good<-complete.cases(temp)
        good<- good & Pop[,yr] > 10e6
# subset with only the complete data rows
        WorldBankCO2<- temp[good,]
  save(WorldBankCO2, file="WorldBankCO2.rda")
}

}
\keyword{datasets}
