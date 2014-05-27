\name{predict.hsmmspec}
\alias{predict.hsmmspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction for hsmmspec}
\description{
  Predicts the underlying state sequence for an observed sequence \code{newdata} given a \code{hsmm} model
}
\usage{
\method{predict}{hsmmspec}(object, newdata, method = "viterbi",M=NA, ...)
}

\arguments{
  \item{object}{ An object of type \code{hsmmspec}}
  \item{newdata}{ A vector or dataframe of observations }
  \item{method}{Prediction method (see details)}
  \item{M}{Maximum number of time spent in a state (truncates the waiting distribution)}   
  \item{\dots}{further arguments passed to or from other methods.}

}
\details{
    If \code{method="viterbi"}, this technique applies the Viterbi
    algorithm for HSMMs, producing the most likely sequence of states
    given the observed data.  If \code{method="smoothed"}, then the
    individually most likely (or smoothed) state sequence is produced,
    along with a matrix with the respective probabilities for each
    state.  This method is different to predict.hsmm in that it takes
    the output from \code{hsmmspec} as input ie. it is useful for people
    who already know their model parameters.
}
\value{
  Returns a \code{hsmm.data} object, suitable for plotting.
  \item{newdata }{A vector or data.frame of observations}
  \item{s }{A vector containing the reconstructed state sequence}
  \item{N }{The lengths of each sequence}
  \item{p }{A matrix where the rows represent time steps and the columns are the probability for the respective state (only produced when \code{method="smoothed"})}
}
\references{ Guedon, Y. (2003),  Estimating hidden semi-Markov chains from discrete sequences, 
	Journal of Computational and Graphical Statistics, Volume 12, Number 3, page 604-639 - 2003}
\author{Jared O'Connell jaredoconnell@gmail.com}

\seealso{ \code{\link{hsmmspec}}, \code{\link{predict.hsmm}}}
\examples{
J <- 3
init <- c(0,0,1)
P <- matrix(c(0,.1,.4,.5,0,.6,.5,.9,0),nrow=J)
B <- list(mu=c(10,15,20),sigma=c(2,1,1.5))
d <- list(lambda=c(10,30,60),shift=c(10,100,30),type='poisson')
model <- hsmmspec(init,P,parms.emission=B,sojourn=d,dens.emission=dnorm.hsmm)
train <- simulate(model,r=rnorm.hsmm,nsim=100,seed=123456)
mean(predict(model,train,M=500)$s!=train$s) #error rate when true model is known
}