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

\name{sim.spatialProcess}
\alias{sim.Krig}
\alias{simSpatialData}
\alias{sim.spatialProcess}
\alias{sim.Krig.approx}
\alias{sim.mKrig.approx}
\alias{simLocal.mKrig}

\title{Conditional simulation of a spatial process}
\description{
Generates exact (or approximate) random draws from the conditional 
distribution of a spatial process given specific observations. This is a 
useful way to characterize the uncertainty in the predicted process from 
data. This is known as conditional simulation in geostatistics or 
generating an ensemble prediction in the geosciences. sim.Krig.grid can 
generate a conditional sample for a large regular grid but is restricted 
to stationary correlation functions. 
 }
\usage{
sim.spatialProcess(object, xp,  M = 1, verbose = FALSE, ...)
simSpatialData(object,  M = 1, verbose = FALSE)

simLocal.mKrig(mKrigObject,  
    predictionGridList = NULL,
    simulationGridList = NULL, 
        gridRefinement = 1, 
                    np = 2, 
                     M = 1,
                    nx = 80,
                    ny = 80, 
               verbose = FALSE,
               delta = NULL, giveWarnings = TRUE,
                          ...)
sim.mKrig.approx(mKrigObject, predictionPoints = NULL,
                 predictionPointsList = NULL, simulationGridList =
                 NULL, gridRefinement = 5, gridExpansion = 1 + 1e-07, M
                 = 1, nx = 40, ny = 40, nxSimulation = NULL,
                 nySimulation = NULL, delta = NULL, verbose = FALSE,...)
 
sim.Krig(object, xp, M = 1, verbose = FALSE, ...)

sim.Krig.approx(object, grid.list = NULL, M = 1, nx = 40, ny = 40,
                 verbose = FALSE, extrap = FALSE,...)


                              

} 
%- maybe also 'usage' for other objects documented here. 
\arguments{
 
 \item{delta}{If the covariance has compact support the simulation method can 
take advantage of this. This is the amount of buffer added for the simulation domain in the circulant embedding method. 
A minimum size would be \code{aRange} for the Wendland but a multiple of this maybe needed to obtain a positive definite 
circulant covariance function.  }

  \item{extrap}{ If FALSE conditional process is not evaluated outside 
     the convex hull of observations. }


  \item{grid.list}{Grid information for evaluating the conditional 
         surface as a grid.list.}
   \item{gridRefinement}{Amount to increase the number of grid points
  for the simulation grid.}

  \item{gridExpansion}{Amount to increase the size of teh simulation
grid. This is used to increase the simulation domain so that the
circulant embedding algorithm works.}

\item{giveWarnings}{If true will warn when more than one observation is in a grid box. This is instead of giving an error and stopping.}

  \item{mKrigObject}{An mKrig Object (or spatialProcess object)}

 
  \item{M}{Number of draws from conditional distribution.}
   \item{np}{Degree of nearest neighbors to use. Default \code{np=2} uses 16
   points in a 4X4 grid for prediction of the off grid point. }
  \item{nx}{ Number of grid points in prediction locations for x coordinate.}
  \item{ny}{ Number of grid points in  prediction locations for x coordinate.}
  

  \item{nxSimulation}{ Number of grid points in the circulant embedding simulation x coordinate.}
  \item{nySimulation}{ Number of grid points in the circulant embedding simulation x coordinate.}

  \item{object}{A Krig object.}
\item{predictionGridList}{A grid list specifying the grid locations for the conditional samples. }
  \item{predictionPoints}{A matrix of locations defining the 
          points for evaluating the predictions.}

   \item{predictionPointsList}{ A \code{grid.list} defining the 
rectangular grid for evaluating the predictions.}

  \item{simulationGridList}{ A \code{gridlist} describing grid for
simulation. If missing this is created from the range of the
locations, \code{nx}, \code{ny}, \code{gridRefinement}, and \code{gridExpansion}
or from the range and and \code{nxSimulation}, \code{nySimulation}.}

 
  \item{xp}{Same as predictionPoints above.}
%  \item{Zp}{The covariate vector or matrix for predicting at the locations xp}

 \item{\dots}{Any other arguments to be passed to the predict function.
 	Usually this is the \code{Z} or \code{drop.Z} argument when there are 
 	additional covariates in the fixed part of the model.
 	(See example below.) }

  \item{verbose}{If true prints out intermediate information. }
}

\details{

These functions generate samples from an unconditional or conditional multivariate (spatial) 
distribution, or an approximate one. The  \strong{unconditional} simulation function,
\code{simSpatialData},
is a handy way to generate synthetic observations from a fitted model.
Typically one would use these for a parametric bootstrap.
The functions that simulate \strong{conditional} distributions are much more involved
in their coding. They are useful
for describing the uncertainty in predictions using the estimated spatial 
process under Gaussian assumptions. An important assumption throughout 
these functions is that all covariance parameters are fixed at their 
estimated or prescribed values from the passed object. Although these functions might
be coded up easily by the users these
versions have the advantage that they take the \code{mKrig}, \code{spatialProcess} or
\code{Krig} objects as a way to specify the
model in an unambiguous way. 

Given a spatial process  h(x)= P(x) + g(x) observed at 

Y.k =  Z(x.k)d + P(x.k) + g(x.k) + e.k

where P(x) is a low order, fixed polynomial and g(x) a Gaussian spatial 
process and Z(x.k) is a vector of covariates that are also indexed by space (such as elevation).
Z(x.k)d is a linear combination of the the covariates with the parameter vector d being a
component of the fixed part of the 
model and estimated in the usual way by generalized least squares.

With Y= Y.1, ..., Y.N,
the goal is to sample the conditional distribution of the process. 
 
[h(x) | Y ]  or the full prediction Z(x)d + h(x)

For fixed a covariance this is just a multivariate normal sampling
problem.  \code{sim.Krig.standard} samples this conditional process at
the points \code{xp} and is exact for fixed covariance parameters.
\code{sim.Krig.grid} also assumes fixed covariance parameters and does
approximate sampling on a grid.

The outline of the algorithm is

0) Find the spatial prediction at the unobserved locations based on
 the actual data. Call this h.hat(x) and this is the conditional mean. 

1) Generate an unconditional spatial process and from this process
simluate synthetic observations. 

2) Use the spatial prediction model ( using the true covariance) to
estimate the spatial process at unobserved locations.

3) Find the difference between the simulated process and its
prediction based on synthetic observations. Call this e(x).

4) h.hat(x) + e(x) is a draw from [h(x) | Y ].

The approximations come int step 1).
Here the field at the observation locations is
approximated using interpolation from the nearest grid points.


\code{sim.spatialProcess} Follows this algorithm exactly. For the case of an
addtional covariate  this of course needs to be included. For a model
with covariates use \code{drop.Z=TRUE} for the function to ignore prediction
using the covariate and generate conditional samples for just the spatial
process and any low order polynomial. Finally, it should be noted that this
function will also work with an \code{mKrig} object because the essential
prediction information in the mKrig and spatialProcess objects are the same.
The naming is through convenience. 

\code{sim.Krig}  Also follows this algorithm exactly but for the older \code{Krig} object.  Note the inclusion of
drop.Z=TRUE or FALSE will determine whether the conditional simulation 
includes the covariates Z or not. (See example below.) 


\code{simNN.mKrig}

\code{sim.Krig.approx} and \code{sim.mKrig.approx} evaluate the
conditional surface on grid and simulates the values of h(x) off the
grid using bilinear interpolation of the four nearest grid
points. Because of this approximation it is important to choose the
grid to be fine relative to the spacing of the observations. The
advantage of this approximation is that one can consider conditional
simulation for large grids -- beyond the size possible with exact
methods. Here the method for simulation is circulant embedding and so
is restricted to stationary fields. The circulant embedding method is
known to fail if the domain is small relative to the correlation
range. The argument \code{gridExpansion} can be used to increase the
size of the domain to make the algorithm work.



}
\value{
\code{sim.Krig and sim.spatialProcess}
 a matrix with rows indexed by the locations in \code{xp} and columns being the 
 \code{M} independent draws.

\code{sim.Krig.approx} a list with components \code{x}, \code{y} and \code{z}. 
x and y define the grid for the simulated field and z is a three dimensional array
 with dimensions \code{c(nx, ny, M)} where the 
first two dimensions index the field and the last dimension indexes the draws.

\code{sim.mKrig.approx} a list with \code{predictionPoints} being the
locations where the field has been simulated.If these have been created from a 
grid list that information is stored in the \code{attributes} of \code{predictionPoints}. 
 \code{Ensemble} is a
matrix where rows index the simulated values of the field and columns
are the different draws, \code{call} is the calling sequence. Not that if \code{predictionPoints}
has been omitted in the call or is created beforehand using \code{make.surface.grid} it is 
easy to reformat the results into an image format for ploting using \code{as.surface}.
e.g. if \code{simOut} is the output object then to plot the 3rd draw:

\preformatted{
     imageObject<- as.surface(simOut$PredictionGrid, simOut$Ensemble[,3] )
     image.plot( imageObject)
}


}

\author{Doug Nychka}
\seealso{ sim.rf, Krig, spatialProcess}
\examples{
## 10 member ensemble for the O3 data

\dontrun{
data( "ozone2")
mKrigObject<- spatialProcess( ozone2$lon.lat, ozone2$y[16,],
                              smoothness=.5)

nx<- 65
ny<- 65

xGridList<- fields.x.to.grid( mKrigObject$x, nx=nx, ny=ny)
xGrid<- make.surface.grid( xGridList)

allTime0<- system.time(
 look0<- sim.spatialProcess(mKrigObject, xp= xGrid, M=10)
)
print( allTime0)

## Local simulation with extra refinement of the grid for embedding
## and same grid size for prediction 
## this runs much faster compared to exact method above 
## as nx, ny are increased  e.g. nx= 128, ny=128 is dramatic difference 

allTime<- system.time(
  look<- simLocal.mKrig(mKrigObject, M=10,nx=nx, ny=ny,
                     gridRefinement = 3,
                     np=3)
) 
print( allTime)
print( look$timing)
}

\dontrun{
## A simple example for setting up a bootstrap 
## M below should be
## set to much larger sample size ( e.g. M <- 200) for better
## statistics

data( ozone2)
obj<- spatialProcess( ozone2$lon.lat,ozone2$y[16,] )
aHat<- obj$summary["aRange"]
lambdaHat<- obj$summary["lambda"]

######### boot strap 
# create M independent copies of the observation vector
# here we just grab the model information from the 
# spatialProcess object above.
#  
# However, one could just create the list 
#  obj<- list( x= ozone2$lon.lat,
#       cov.function.name="stationary.cov",
#     summary= c( tau= 9.47, sigma2= 499.79, aRange= .700),
#    cov.args= list( Covariance="Matern", smoothness=1.0),
#     weights= rep( 1, nrow(ozone2$lon.lat) )
# )
# Here summary component  has the parameters
# tau, sigma2 and aRange
# and cov.args component has the remaining ones.

set.seed(223)
M<- 25
ySynthetic<- simSpatialData( obj, M)

bootSummary<- NULL

for(  k in 1:M){
cat( k, " ")
# here the MLEs are found using the easy top level level wrapper
# see mKrigMLEJoint for a more efficient strategy
  newSummary<- spatialProcess(obj$x,ySynthetic[,k],
                    cov.params.start= list(
			                   aRange = aHat,
			                  lambda = lambdaHat)
                               )$summary
  bootSummary<- rbind( bootSummary, newSummary)
  }
cat( fill= TRUE)
# the results and 95% confidence interval  
  stats( bootSummary )

  obj$summary
  tmpBoot<- bootSummary[,c("lambda", "aRange") ]
  confidenceInterval <- apply(tmpBoot, 2,
                               quantile, probs=c(0.025,0.975) )
# compare to estimates used as the "true" parameters			       
  obj$summary[2:5] 
  print( t(confidenceInterval) )
# compare to confidence interval using large sample theory  
  print( obj$CITable)
}

\dontrun{
# conditional simulation with covariates
# colorado climate example
  data(COmonthlyMet)
  fit1E<- spatialProcess(CO.loc,CO.tmin.MAM.climate, Z=CO.elev   )
# conditional simulation at missing data
  good<- !is.na(CO.tmin.MAM.climate ) 
  infill<- sim.spatialProcess( fit1E, xp=CO.loc[!good,], 
                Z= CO.elev[!good], M= 10)
# get an elevation grid  ... NGRID<- 50 gives a nicer image but takes longer 
 NGRID <- 25  
 # get elevations on a grid  
   COGrid<- list( x=seq( -109.5, -101, ,NGRID), y= seq(39, 41.5,,NGRID) )
   COGridPoints<- make.surface.grid( COGrid)
 # elevations are a bilinear interpolation from the 4km
 # Rocky Mountain elevation fields data set.   
   data( RMelevation)
   COElevGrid<- interp.surface( RMelevation, COGridPoints )
# NOTE call to sim.spatialProcess treats the grid points as just a matrix
# of locations the plot has to "reshape" these into a grid 
# to use with image.plot 
   SEout<- sim.spatialProcess( fit1E, xp=COGridPoints,  Z= COElevGrid, M= 30)
# for just the smooth surface in lon/lat
#  SEout<- sim.spatialProcess( fit1E, xp=COGridPoints,  drop.Z=TRUE, M= 30)
# in practice M should be larger to reduce Monte Carlo error.      
   surSE<- apply( SEout, 2, sd )
   image.plot( as.surface( COGridPoints, surSE)) 
   points( fit1E$x, col="magenta", pch=16) 
   
}

data( ozone2)
set.seed( 399)
# fit to day 16 from Midwest ozone data set.
  out<- Krig( ozone2$lon.lat, ozone2$y[16,], Covariance="Matern", 
            aRange=1.0,smoothness=1.0, na.rm=TRUE)

# NOTE aRange =1.0 is not the best choice but 
# allows the sim.rf circulant embedding algorithm to 
# work without increasing the domain.

#six missing data locations
 xp<-  ozone2$lon.lat[ is.na(ozone2$y[16,]),]

# 5 draws from process at xp given the data 
# this is an exact calculation
 sim.Krig( out,xp, M=5)-> sim.out

# Compare: stats(sim.out)[3,] to  Exact: predictSE( out, xp)
# simulations on a grid
# NOTE this is approximate due to the bilinear interpolation
# for simulating the unconditional random field. 
# also more  grids points ( nx and  ny) should be used  

sim.Krig.approx(out,M=5, nx=20,ny=20)-> sim.out

# take a look at the ensemble members. 

predictSurface( out, grid= list( x=sim.out$x, y=sim.out$y))-> look

zr<- c( 40, 200)

set.panel( 3,2)
image.plot( look, zlim=zr)
title("mean surface")
for ( k in 1:5){
image( sim.out$x, sim.out$y, sim.out$z[,,k], col=tim.colors(), zlim =zr)
}



\dontrun{
data( ozone2)
y<- ozone2$y[16,]
good<- !is.na( y)
y<-y[good]
x<- ozone2$lon.lat[good,]
O3.fit<- mKrig( x,y, Covariance="Matern", aRange=.5,smoothness=1.0, lambda= .01 )
set.seed(122)
O3.sim<- sim.mKrig.approx( O3.fit, nx=100, ny=100, gridRefinement=3, M=5 )
set.panel(3,2)
surface( O3.fit)
for ( k in 1:5){
image.plot( as.surface( O3.sim$predictionPoints, O3.sim$Ensemble[,k]) )
}
# conditional simulation at missing data
xMissing<- ozone2$lon.lat[!good,]
O3.sim2<- sim.mKrig.approx( O3.fit, xMissing, nx=80, ny=80,
                            gridRefinement=3, M=4 )
}
\dontrun{
#An example for fastTps:
  data(ozone2)
  y<- ozone2$y[16,]
  good<- !is.na( y)
  y<-y[good]
  x<- ozone2$lon.lat[good,]
  O3Obj<- fastTps( x,y, aRange=1.5 )
# creating a quick grid list based on ranges of locations
  grid.list<- fields.x.to.grid( O3Obj$x, nx=100, ny=100)
# controlling the grids
  xR<- range( x[,1], na.rm=TRUE)
  yR<- range( x[,2], na.rm=TRUE)
  simulationGridList<- list( x= seq(xR[1],xR[2],,400),
         y= seq( yR[1],yR[2], ,400))
# very fine localized prediction grid
    O3GridList<- list( x= seq( -90.5,-88.5,,200), y= seq( 38,40,,200))
    O3Sim<- sim.mKrig.approx( O3Obj, M=5, predictionPointsList=O3GridList,
                  simulationGridList = simulationGridList)
}
}
\keyword{spatial}
% at least one, from doc/KEYWORDS
