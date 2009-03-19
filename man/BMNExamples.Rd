\name{BMNExamples}
\alias{BMNSamples}
\title{Sampling data using Gibbs sampling for use in the examples}
\description{
    Sampling from pairwise binary Markov model using Gibbs sampling. This function is not efficient and only intended to be used in the examples.
}
\usage{
BMNSamples(Theta, numSamples, burnIn, skip)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Theta}{Parameter matrix for the model from which the data is being generated.}
  \item{numSamples}{Number of samples to return.}
  \item{burnIn}{Number of samples to discard as burn in.}
  \item{skip}{Number of samples to discard in-between returned samples.}
}
\details{
\code{BMNSamples} generates \code{numSamples} by using Gibbs sampling. When using Gibbs sampling, it is necessary to discard the initial samples, which is controlled by the parameter \code{burnIn}. In order for the drawn samples to be independent, samples in-between also have to be discarded, which is controlled by \code{skip}.
}
\value{
Returns a matrix of 0 and 1 of size \code{numSamples} times \code{p} where \code{p} is the number of rows of \code{Theta}. 
}
\seealso{\code{\link{BMNPseudo}}, \code{\link{BMNExact}}}
\author{Holger Hoefling}
\keyword{multivariate}
