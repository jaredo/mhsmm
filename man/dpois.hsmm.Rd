\name{dpois.hsmm}
\alias{dpois.hsmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Emission density function for Poisson emission distribution }
\description{
  Calculates the density of observations \code{x} for state \code{j} given the parameters in \code{model}.  This is used for
  a Poisson emission distribution of a HMM or HSMM and is a suitable prototype for user's to make their own custom distributions.
}
\usage{
dpois.hsmm(x, j, model)
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

\seealso{ \code{\link{mstep.pois}}, \code{\link{rpois.hsmm}}}

\author{Jared O'Connell jaredoconnell@gmail.com}
\examples{
  J<-3
  initial <- rep(1/J,J)
  P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
  b <- list(lambda=c(1,3,6))
  model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dpois.hsmm)
  model
  train <- simulate(model, nsim=300, seed=1234, rand.emis=rpois.hsmm)
  plot(train,xlim=c(0,100))  
  h1 = hmmfit(train,model,mstep=mstep.pois)
}