\name{offGridWeights}
\alias{offGridWeights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weights to predict off grid locations from a rectangular  grid using nearest neighbors and Kriging. 
%%  ~~function to do ... ~~
}
\description{
Based on a stationary Gaussian process model creates a sparse matrix to predict off grid values (aka interpoltate) from an equally spaced rectangular grid. The sparsity comes about because only a fixed number of neighboring grid points (np) are used in the prediction.  The prediction variance is also give in the returned object. This function is used as the basis for approximate conditional simulation for large spatial datasets. 

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
offGridWeights(s, gridList, np = 2, mKrigObject = NULL, Covariance = NULL, covArgs = NULL, aRange = NULL, sigma2 = NULL,giveWarnings = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{s}{
  Off grid spatial locations
%%     ~~Describe \code{s} here~~
}
  \item{gridList}{
  A list as the gridList format ( x and y components) that describes the rectagular grid. The grid must have at least np extra grid points beyond the range of the points in \code{s}
%%     ~~Describe \code{gridList} here~~
}
  \item{np}{
  Number of nearest neighbor grid points to use for prediction.  \code{np = 1} will use the 4 grid points that bound the off grid point.  \code{np = 2} will be a 4X4 subgrid with the middle grid box containing the off gird point. In general there will be  \code{(2*np)^2} neighboring points uses. 
%%     ~~Describe \code{np} here~~
}
  \item{mKrigObject}{
  The output object (aka list) from either mKrig or spatialProcess. This has the information about the covariance function used to do the Kriging. 
  The following items are coded in place of not supplying this object.  See the example below for more details.
%%     ~~Describe \code{mKrigObject} here~~
}
  \item{Covariance}{
  The stationary covariance function  (taking pairwise distances as its first argument.)
%%     ~~Describe \code{Covariance} here~~
}
  \item{covArgs}{
  If \code{mKrigObject} is not specified  a list giving any additional arguments for the covariance function. 
%%     ~~Describe \code{covArgs} here~~
}
  \item{aRange}{
  The range parameter. 
  
%%     ~~Describe \code{aRange} here~~
}
  \item{sigma2}{
  Marginal variance of the process. 
%%     ~~Describe \code{sigma2} here~~
}
 \item{giveWarnigs}{If TRUE will warn if two or more observations
 are in the same grid box.
}
}
\details{
This function creates the interpolation weights taking advantage of some
efficiency in the covariance function being stationary, use
of a fixed configuration of nearest neighbors, and Kriging
predictions from a rectangular grid. 

The returned matrix is in spam sparse matrix format. See
example below for the "one-liner" to make the prediction
once the weights are computed. Although created primarily
for conditional simulation of a spatial process this
function is also useful for interpolating to off grid
locations from a rectangular field. 

The user is also referred to the testing script
\code{offGridWeights.test.R} in \code{tests} where the 
Kriging is computed explicitly and tested against the sparse
matrix computation. This is helpful in defining exactly what
is being computed.  

%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{B}{A sparse matrix that is of dimension mXn with m the number of locations (rows) in \code{s} and n being the total number of grid points. \code{n = length(gridList$x)*length(gridList$y) }
}
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{predictionVariance}{A vector of length as the rows of \code{s} with the Kriging prediction variance based on the nearest neighbor prediction and the specified covariance function. }
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{ Douglas Nychka and Maggie Bailey
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{interp.surface} 
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# an M by M  grid
M<- 400
xGrid<- seq( -1, 1, length.out=M)
gridList<- list( x= xGrid,
                 y= xGrid
                 )
 np<- 3 
 n<- 100
# sample n locations but avoid margins 
set.seed(123)
s<- matrix(    runif(n*2, xGrid[(np+1)],xGrid[(M-np)]),
               n, 2 )

                  
obj<- offGridWeights( s, gridList, np=3,
                   Covariance="Matern",
                   aRange = .1, sigma2= 1.0,
                   covArgs= list( smoothness=1.0)
                   )
# make the predictions  by obj$B%*%c(y)
# where y is the matrix of values on the grid
 
# try it out on a simulated  Matern field  
CEobj<- circulantEmbeddingSetup( gridList,  
                  cov.args=list(
                  Covariance="Matern",
                   aRange = .1,
                    smoothness=1.0)
                    )
 set.seed( 333)                   
Z<- circulantEmbedding(CEobj)

#
# Note that grid values are "unrolled" as a vector
# for multiplication
# predOffGrid<- obj$B%*% c( Z)

predOffGrid<- obj$B\%*\% c( Z)

set.panel( 1,2)
zr<- range( c(Z))
image.plot(gridList$x, gridList$y, Z, zlim=zr)
bubblePlot( s[,1],s[,2], z= predOffGrid , size=.5,
highlight=FALSE, zlim=zr)
 set.panel()
 
                  
  
}
\keyword{spatial}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
