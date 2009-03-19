\name{BMNJT}
\alias{BMNJT}
\title{ Junction tree algorithm for pairwise binary Markov networks }
\description{
Calculates the expections and second moments for a pairwise binary Markov network using the junction tree algorithm.
}
\usage{
BMNJT(thetaMat, adjMat = NULL, var = NULL, onlyActive = FALSE, timeout = 60)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{thetaMat}{Theta matrix of the model (assumed to be symmetric).}
  \item{adjMat}{Adjacency matrix of the underlying graph. If not given, inferred from \code{thetaMat}. If given, values in \code{thetaMat} that do not correspond to an edge in the graph are being ignored.}
  \item{var}{if given, the expectation and second moments of variable var only will be calculated. Otherwise, for all variables.}
  \item{onlyActive}{If \code{TRUE}, only values in the second moment matrix corresponding to edges in the graph will be calculated (much faster). Otherwise, all entries will be calculated.}
  \item{timeout}{Number of seconds until the functions terminates automatically.}
}
\details{
Uses the junction tree algorithm to calculate the matrix of second moments of the underlying pairwise binary Markov model. This is useful for inference on this type of models as the matrix of second moments is the derivative of the partition function in this class of models. The function provides the option to calculate the second moments only for one variable with all the others (variables are numbered from 1 to p; see \code{var}). Also, the juntion tree algorithm is by a factor p more efficient if only second moments corresponding to edges in the underlying graph (non-zero values in the \code{thetaMat} or \code{adjMat}) are being calculated. This can be done using the switch \code{onlyActive}.
}
\value{
The return value is a list with elements:
  \item{Expectation}{Expectation of the variables.}
  \item{SecondMomentMatrix}{Second moments of the variables (only present if \code{Var=NULL}).}
  \item{SecondMomentVector}{Vector of second moments of variable \code{Var} with the others if \code{Var != NULL}.}
}
\seealso{\code{\link{BMNPseudo}}, \code{\link{BMNExact}}}
\examples{
library(BMN)
Theta = matrix(numeric(25), ncol=5)
Theta[1,1]=0.5; Theta[2,2]=0.5; Theta[3,3]=0; Theta[4,4]= -0.5; Theta[5,5]= 0.5;
Theta[1,2]=Theta[2,1]=1; Theta[1,4]=Theta[4,1]=1; Theta[2,3]=Theta[3,2]= -1;

BMNJT(Theta)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\author{Holger Hoefling}
\keyword{multivariate}
