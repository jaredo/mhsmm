\name{simulate.hsmmspec}
\alias{simulate.hsmmspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulation for HSMMs }
\description{
  Simulates values for a specified hidden semi-Markov model
}
\usage{
\method{simulate}{hsmmspec}(object, nsim, seed = NULL,rand.emission=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{hsmmspec} object }
  \item{nsim}{An integer or vector of integers (for multiple sequences) specifying the length of the sequence(s) }
  \item{seed}{ \code{seed} for the random number generator}
  \item{rand.emission}{The function used to generate observations from the emission distribution}
  \item{\dots}{further arguments passed to or from other methods.}

}
\details{
  If \code{nsim} is a single integer then a HSMM of that length is produced.  If \code{nsim} is a vector of integers, then \code{length(nsim)} sequences are generated with respective lengths.
}
\value{
  An object of class hmmdata
  \item{x}{A vector of length \code{sum(N)} - the sequence(s) of observed values}
  \item{s}{A vector of length \code{sum(N)} - the sequence(s) of hidden states}
  \item{N}{A vector of the length of each observation sequence (used to segment x and s)}
}

\references{ Guedon, Y. (2003),  Estimating hidden semi-Markov chains from discrete sequences, 
	Journal of Computational and Graphical Statistics, Volume 12, Number 3, page 604-639 - 2003 }
\author{ Jared O'Connell jaredoconnell@gmail.com }


\seealso{ 
  \code{\link{hsmmfit}},
  \code{\link{hsmmspec}},
  \code{\link{predict.hsmm}}
}

\examples{
J <- 3
init <- c(0,0,1)
P <- matrix(c(0,.1,.4,.5,0,.6,.5,.9,0),nrow=J)
B <- list(mu=c(10,15,20),sigma=c(2,1,1.5))
d <- list(lambda=c(10,30,60),shift=c(10,100,30),type='poisson')
model <- hsmmspec(init,P,parms.emission=B,sojourn=d,dens.emission=dnorm.hsmm)
train <- simulate(model,rand.emis=rnorm.hsmm,nsim=100,seed=123456)
plot(train,xlim=c(0,400))
}