\name{BMNExact}
\alias{BMNExact}
\alias{BMNExact.single}
\title{Exact inference in L1-penalized Binary Markov Model}
\description{
    Functions to calculate a parameter estimates for L1-penalized Binary Markov Models.
}
\usage{
BMNExact(X, rhoVec, thrCheck=1e-3, thrPseudo=1e-5, ThetaStart=NULL, 
verbose = FALSE, maxIter=100, timeout=60, penalize.diag=FALSE)
BMNExact.single(X, rho, thrCheck=1e-3, thrPseudo=1e-5, ThetaStart=NULL, 
verbose = FALSE, maxIter=100, timeout=60, penalize.diag=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Input data matrix consisting of 0-1 entries. Has n rows and p columns.}
  \item{rho}{Value of the penalty parameter; If a non-negative p-by-p matrix is given, it is used as the penalty structure.}
  \item{rhoVec}{Gives all values of rho for which the solution should be calculated.}
  \item{thrCheck}{Error threshold at which convergence is declared.}
  \item{thrPseudo}{Error threshold for the interal pseudolikelihood algorithm.}
  \item{ThetaStart}{Starting value for Theta, has to be a p-by-p matrix.}
  \item{verbose}{Should status messages be printed.}
  \item{maxIter}{Maximum number of iteratios to run.}
  \item{timeout}{Number of seconds after which the procedure is stopped; for the path algorithm, this is reset for every value of rho.}
  \item{penalize.diag}{Should the diagonal be penalized?}
}
\details{
The function \code{BMNExact} fits a penalized pairwise binary Markov model to the data provided as matrix \code{X} for each of the elements in the penalty parameter vector \code{rhoVec} (note that \code{rhoVec} will be sorted in increasing order). Internally, the function \code{BMNExact.single} is called for each entry in \code{rhoVec} and the results are collected as described below. 
}
\value{
\item{rho}{Vector of penalty parameters sorted in increasing order.}
\item{ThetaList}{A list of Theta pxp matrices, corresponding to the penalty parameters in rho.}
\item{success}{A logical vector of the same length as rho. True, if the function succeeded for the corresponding value in rho.}
\item{penalize.diag}{Logical. Indicates if the diagonal was penalized (same as input value \code{penalize.diag}.}
}
\seealso{\code{\link{BMNPseudo}}, \code{\link{BMNJT}}}
\examples{
library(BMN)
Theta = matrix(numeric(25), ncol=5);
Theta[1,1]=0.5; Theta[2,2]=0.5; Theta[3,3]=0; Theta[4,4]= -0.5; Theta[5,5]= 0.5;
Theta[1,2]=Theta[2,1]=1; Theta[1,4]=Theta[4,1]=1; Theta[2,3]=Theta[3,2]= -1;
numSamples=1000; burnIn=100; skip=1;

simData = BMNSamples(Theta, numSamples, burnIn, skip)

rhoVec = c(0.01, 0.02, 0.03)

exactPath = BMNExact(simData, rhoVec)
exactSingle = BMNExact.single(simData, 0.02)
}
\author{Holger Hoefling}
\keyword{multivariate}
