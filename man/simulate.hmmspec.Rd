\name{simulate.hmmspec}
\alias{simulate.hmmspec}
\title{Simulation of hidden Markov models}
\description{
Simulates data from a hidden Markov model
}
\usage{
\method{simulate}{hmmspec}(object, nsim, seed = NULL, rand.emission=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{hmmspec} object }
  \item{nsim}{An integer or vector of integers (for multiple sequences)
    specifying the length of the sequence(s) } 
  \item{seed}{ \code{seed} for the random number generator}
  \item{rand.emission}{The function used to generate observations from the emission distribution}
  \item{\dots}{further arguments passed to or from other methods.}
  
}
\details{
  If \code{nsim} is a single integer then a HMM of that length is
  produced.  If \code{nsim} is a vector of integers, then
  \code{length(nsim)} sequences are generated with respective lengths. 
}
\value{
  An object of class hmmdata
  \item{x}{A vector of length \code{sum(N)} - the sequence(s) of observed values}
  \item{s}{A vector of length \code{sum(N)} - the sequence(s) of hidden states}
  \item{N}{A vector of the length of each observation sequence (used to segment x and s)}
}
\references{Rabiner, L. (1989), A tutorial on hidden Markov models and selected applications
in speech recognition, Proceedings of the IEEE, 77, 257-286.}
\author{ Jared O'Connell jaredoconnell@gmail.com }

\seealso{ \code{\link{hmmspec}}, \code{link{predict.hmm}} }

\examples{

J<-3
initial <- rep(1/J,J)
P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
b <- list(mu=c(-3,0,2),sigma=c(2,1,.5))
model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dnorm.hsmm)
train <- simulate(model, nsim=100, seed=1234, rand.emis=rnorm.hsmm)
plot(train)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
