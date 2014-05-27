
\name{hmmspec}
\alias{hmmspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Specificatin of HMMs }
\description{
  Creates a model specficiation for a hidden Markov model
}
\usage{
hmmspec(init, trans, parms.emission, dens.emission, rand.emission=NULL,mstep=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{Distribution of states at t=1 ie. P(S=s) at t=1}
  \item{trans}{The transition matrix of the Markov chain}
  \item{parms.emission}{A list containing the parameters of the emission distribution}
  \item{dens.emission}{Density function of the emission distribution.}
  \item{rand.emission}{The function used to generate observations from the emission distribution}
  \item{mstep}{ Re-estimates the parameters of density function on each iteration}
  
}
\value{
  A hmmspec object
}
\references{Jared O'Connell, Soren Hojsgaard (2011). Hidden Semi Markov Models for
Multiple Observation Sequences: The mhsmm Package for R., 
Journal of Statistical Software, 39(4), 1-22., URL
http://www.jstatsoft.org/v39/i04/.

Rabiner, L. (1989), A tutorial on hidden Markov models and selected applications
in speech recognition, Proceedings of the IEEE, 77, 257-286.}
\author{ Jared O'Connell jaredoconnell@gmail.com }

\seealso{ \code{\link{simulate.hmmspec}}, \code{\link{simulate.hmmspec}}, \code{\link{hmmfit}}, \code{\link{predict.hmm}}}

