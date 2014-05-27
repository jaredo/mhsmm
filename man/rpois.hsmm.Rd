\name{rpois.hsmm}
\alias{rpois.hsmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Random number generation from a Poisson distributed emission distribution }
\description{
  This generates values from a Poisson distributed emission state \code{j} given parameters in \code{model}.
}
\usage{
rpois.hsmm(j, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{ An integer representing the state }
  \item{model}{ A \code{hmmspec} or \code{hsmmspec} object }
}
\details{
  This is essentially a wrapper for \code{rpois}.  Users may build functions with the same 
  arguments and return values so they can use their own custom emission distributions.
}
\value{
  A single value from the emission distribution.
}

\seealso{ \code{\link{mstep.pois}}, \code{\link{dpois.hsmm}}}

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