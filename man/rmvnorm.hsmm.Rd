\name{rmvnorm.hsmm}
\alias{rmvnorm.hsmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Random number generation from a multivariate normal distributed emission distribution }
\description{
  This generates values from a multivariate normal distributed emission
  state \code{j} given parameters in \code{model}. 
}
\usage{
rmvnorm.hsmm(j, model)
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
\seealso{ \code{\link{dmvnorm.hsmm}}, \code{\link{mstep.mvnorm}} }

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

