\name{mstep.mvnorm}
\alias{mstep.mvnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performs re-estimation (the M-step) for a multivariate normal emission distribution}
\description{
  Re-estimates the parameters of a multivariate normal emission
  distribution as part of the EM algorithm for HMMs and HSMMs.  
  This is called by the \code{hmm} and \code{hsmm} functions.  It is a
  suitable prototype function for users wishing to design their own
  emission distributions. 
}
\usage{
mstep.mvnorm(x, wt)
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
  \item{mu}{A list of length J contain the mean vectors}
  \item{sigma}{A list of length J containing the covariance matrices}
}

\seealso{ \code{\link{dmvnorm.hsmm}}, \code{\link{rmvnorm.hsmm}} }

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

