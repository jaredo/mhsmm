\name{dmvnorm.hsmm}
\alias{dmvnorm.hsmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Emission ensity function for a multivariate normal emission distribution }
\description{
  Calculates the density of observations \code{x} for state \code{j} given the parameters in \code{model}.  This is used for
  a multivariate Gaussian emission distribution of a HMM or HSMM and is a suitable prototype for user's to make their own custom distributions.
}
\usage{
dmvnorm.hsmm(x, j, model)
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
\seealso{
  \code{\link{mstep.mvnorm}},
  \code{\link{rmvnorm.hsmm}}
}

\author{Jared O'Connell jaredoconnell@gmail.com}
\examples{
  J<-2
  initial <- rep(1/J,J)
  P <- matrix(c(.3,.5,.7,.5),nrow=J)
  b <- list(mu=list(c(-3,0),c(1,2)),sigma=list(diag(2),matrix(c(4,2,2,3), ncol=2)))
  model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dmvnorm.hsmm)
  model
  train <- simulate(model, nsim=300, seed=1234, rand.emis=rmvnorm.hsmm)
  plot(train,xlim=c(0,100))
  h1 = hmmfit(train,model,mstep=mstep.mvnorm)
}

