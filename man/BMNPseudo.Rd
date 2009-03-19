\name{BMNPseudo}
\alias{BMNPseudo}
\alias{BMNPseudo.single}
\title{Pseudo-likelihood inference in L1-penalized Binary Markov Model}
\description{
    Functions to calculate approximate parameter estimates for L1-penalized Binary Markov Models.
}
\usage{
BMNPseudo(X, rhoVec, Delta=NULL, ThetaStart=NULL, 
maxError=1e-5, verbose=FALSE, maxIter=100, penalize.diag=FALSE, 
stepSize=1, performLineSearch=FALSE)
BMNPseudo.single(X, rho, Delta=NULL, ThetaStart=NULL, 
maxError=1e-5, maxIter=100, penalize.diag=FALSE, 
stepSize=1, performLineSearch=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{Input data matrix consisting of 0-1 entries. Has n rows and p columns.}
  \item{rho}{Value of the penalty parameter; If a non-negative p-by-p matrix is given, it is used as the penalty structure.}
  \item{rhoVec}{Gives all values of rho for which the solution should be calculated.}
  \item{Delta}{Adjustmeant to the gradient.}
  \item{ThetaStart}{Starting value for Theta, has to be a p-by-p matrix.}
  \item{maxError}{convergence threshold for the algorithm.}
  \item{verbose}{Print status messages.}
  \item{maxIter}{Maximum number of iteratios to run.}
  \item{penalize.diag}{Should the diagonal be penalized?}
  \item{stepSize}{Stepsize of the algorithm; should be 1 or less.}
  \item{performLineSearch}{If TRUE, a line search is performed; takes longer but is guarateed to converge.}
}
\details{
The function \code{BMNPseudo} fits an approximate penalized pairwise binary Markov model to the data provided as matrix \code{X} for each of the elements in the penalty parameter vector \code{rhoVec} (note that \code{rhoVec} will be sorted in increasing order). Internally, the function \code{BMNExact.single} is called for each entry in \code{rhoVec} and the results are collected as described below. 
}
\value{
\item{rho}{Vector of non-negative penalty parameters sorted in decreasing order.}
\item{ThetaList}{A list of Theta pxp matrices, corresponding to the penalty parameters in rho.}
\item{success}{A logical vector of the same length as rho. True, if the function succeeded for the corresponding value in rho.}
\item{penalize.diag}{Logical. Indicates if the diagonal was penalized (same as input value \code{penalize.diag}.}
}
\seealso{\code{\link{BMNExact}}, \code{\link{BMNJT}}}
\examples{
library(BMN)
Theta = matrix(numeric(25), ncol=5);
Theta[1,1]=0.5; Theta[2,2]=0.5; Theta[3,3]=0; Theta[4,4]= -0.5; Theta[5,5]= 0.5;
Theta[1,2]=Theta[2,1]=1; Theta[1,4]=Theta[4,1]=1; Theta[2,3]=Theta[3,2]= -1;
numSamples=1000; burnIn=100; skip=1;

simData = BMNSamples(Theta, numSamples, burnIn, skip)

rhoVec = c(0.01, 0.02, 0.03)

pseudoPath = BMNPseudo(simData, rhoVec)
pseudoSingle = BMNPseudo.single(simData, 0.02)
}
\author{Holger Hoefling}
\keyword{multivariate}
