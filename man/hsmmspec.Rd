\name{hsmmspec}
\alias{hsmmspec}
\title{Hidden semi-Markov model specification}
\description{
  Creates a model specification of a hidden semi-Markov model.
}
\usage{
hsmmspec(init,transition,parms.emission,sojourn,dens.emission,
	rand.emission=NULL,mstep=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{init}{Distribution of states at t=1 ie. P(S=s) at t=1}
  \item{transition}{The transition matrix of the embedded Markov chain (diagonal must be 0)}
  \item{parms.emission}{A list containing the parameters of the emission distribution}
  \item{sojourn}{A list containining the parameters and type of sojourn distribtuion (see Details)}
  \item{dens.emission}{ Density function of the emission distribution }
  \item{rand.emission}{The function used to generate observations from the emission distribution}              
  \item{mstep}{ Re-estimates the parameters of density function on each iteration }
}
\details{
  The sojourn argument provides a list containing the parameters for the available sojourn distributions.  Available sojourn distributions are
  shifted Poisson, Gamma and non-parametric.
  
  In the case of the Gamma distribution, sojourn is a list with vectors shape and scale (the Gamma parameters in dgamma), both of length J.  Where J is the number of states. See \code{reprocows} for an example using Gamma sojourn distributions.
  
  In the case of shifted Poisson, sojourn is list with vectors shift and lambda, both of length J.  See \code{hsmmfit} for an example using shifted Poisson sojourn distributions.
  
  In the case of non-parametric, sojourn is a list containing a M x J matrix.  Where entry (i,j) is the probability of a sojourn of length i in state j.  See \code{hsmmfit} for an example using shifted non-parametric sojourn distributions.
  
}
\value{
  An object of class hsmmspec
}
\references{Jared O'Connell, Soren Hojsgaard (2011). Hidden Semi Markov Models for
Multiple Observation Sequences: The mhsmm Package for R., 
Journal of Statistical Software, 39(4), 1-22., URL
http://www.jstatsoft.org/v39/i04/.

Guedon, Y. (2003),  Estimating hidden semi-Markov chains from discrete sequences, 
	Journal of Computational and Graphical Statistics, Volume 12, Number 3, page 604-639 - 2003}
\author{ Jared O'Connell jaredoconnell@gmail.com }
\seealso{ \code{\link{hsmmfit}} ,\code{\link{simulate.hsmmspec}}, \code{\link{predict.hsmm}}}
