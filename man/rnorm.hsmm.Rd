\name{rnorm.hsmm}
\alias{rnorm.hsmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Random number generation from a normally distributed emission distribution }
\description{
  This generates values from a normally distributed emission state \code{j} given parameters in \code{model}.
}
\usage{
rnorm.hsmm(j, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{ An integer representing the state }
  \item{model}{ A \code{hmmspec} or \code{hsmmspec} object }
}
\details{
  This is essentially a wrapper for \code{rnorm}.  Users may build functions with the same 
  arguments and return values so they can use their own custom emission distributions.
}
\value{
  A single value from the emission distribution.
}

\author{Jared O'Connell jaredoconnell@gmail.com}
