\name{dnorm.hsmm}
\alias{dnorm.hsmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Emission density function for normal emission distribution }
\description{
  Calculates the density of observations \code{x} for state \code{j} given the parameters in \code{model}.  This is used for
  the Gaussian emission distribution of a HMM or HSMM and is a suitable prototype for user's to make their own custom distributions.
}
\usage{
dnorm.hsmm(x, j, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Observed value }
  \item{j}{ State }
  \item{model}{ A \code{hsmmspec} or \code{hmmspec} object }
}
\details{
  This is used by \code{hmm} and \code{hsmm} to calculate densities for use in the E-step of the EM algorithm.  
  It can also be used as a template for users wishing to building their own emission distributions
}
\value{
  A vector of probability densities.
}
\author{Jared O'Connell jaredoconnell@gmail.com}
