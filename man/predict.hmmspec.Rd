\name{predict.hmmspec}
\alias{predict.hmmspec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prediction function for hmmspec }
\description{
  Predicts the underlying state sequence for an observed sequence \code{newdata} given a \code{hmmspec} model
}
\usage{
\method{predict}{hmmspec}(object, newdata,method = "viterbi", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class \code{hmm}}
  \item{newdata}{ A vector or list of observations }
  \item{method}{Prediction method (see details)}
  \item{\dots}{further arguments passed to or from other methods.}

}
\details{
    If \code{method="viterbi"}, this technique applies the Viterbi
    algorithm for HMMs, producing the most likely sequence of states
    given the observed data.  If \code{method="smoothed"}, then the
    individually most likely (or smoothed) state sequence is produced,
    along with a matrix with the respective probabilities for each
    state.  This function differs from predict.hmm in that it takes the
    output from hmmspec ie. this is useful when users already know their
    parameters and wish to make predictions.
}
\value{
  Returns a \code{hsmm.data} object, suitable for plotting.
  \item{newdata }{A vector or list of observations}
  \item{s }{A vector containing the reconstructed state sequence}
  \item{N }{The lengths of each sequence}
  \item{p }{A matrix where the rows represent time steps and the columns are the probability for the respective state (only produced when \code{method="smoothed"})}  
}
\references{Rabiner, L. (1989), A tutorial on hidden Markov models and selected applications
in speech recognition, Proceedings of the IEEE, 77, 257-286.}

\author{Jared O'Connell jaredoconnell@gmail.com}
\seealso{\code{hmmspec}}
\examples{
     J<-3
     initial <- rep(1/J,J)
     P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
     b <- list(mu=c(-3,0,2),sigma=c(2,1,.5))
     model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dnorm.hsmm)    
     train <- simulate(model, nsim=300, seed=1234, rand.emis=rnorm.hsmm)
     mean(predict(model,train)$s!=train$s) #error rate when true model is known
}