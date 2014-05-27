\name{mstep.norm}
\alias{mstep.norm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performs re-estimation (the M-step) for a normal emission distribution}
\description{
  Re-estimates the parameters of a normal emission distribution as part of the EM algorithm for HMMs and HSMMs.
  This is called by the \code{hmm} and \code{hsmm} functions.  It is a suitable prototype function for users
  wishing to design their own emission distributions.
}
\usage{
mstep.norm(x, wt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A vector of observed values}
  \item{wt}{ A T x J matrix of weights.  Column entries are the weights for respective states. }
}
\details{
  Users may write functions that take the same arguments and return the same values for their
  own custom emission distributions.
}
\value{
  Returns the \code{emission} slot of a \code{hmmspec} or \code{hsmmspec} object
  \item{mu}{Vector of length J contain the means}
  \item{sigma}{Vector of length J containing the variances}
}

\author{Jared O'Connell jaredoconnell@gmail.com}
