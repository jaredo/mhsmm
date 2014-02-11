\name{predict.hmm}
\alias{predict.hmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Prediction function for hmm }
\description{
  Predicts the underlying state sequence for an observed sequence \code{newdata} given a \code{hmm} model
}
\usage{
\method{predict}{hmm}(object, newdata,method = "viterbi", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of class \code{hmm}}
  \item{newdata}{ A vector or data.frame of observations }
  \item{method}{Prediction method (see details)}
  \item{\dots}{further arguments passed to or from other methods.}

}
\details{
    If \code{method="viterbi"}, this technique applies the Viterbi algorithm for HMMs, producing the most likely sequence of states given the observed data.  If \code{method="smoothed"}, then the individually most likely (or smoothed) state sequence is produced, along with a matrix with the respective probabilities for each state.
}
\value{
  Returns a \code{hsmm.data} object, suitable for plotting.
  \item{newdata }{A vector or data.frame of observations}
  \item{s }{A vector containing the reconstructed state sequence}
  \item{N }{The lengths of each sequence}
  \item{p }{A matrix where the rows represent time steps and the columns are the probability for the respective state (only produced when \code{method="smoothed"})}  
}
\references{Rabiner, L. (1989), A tutorial on hidden Markov models and selected applications
in speech recognition, Proceedings of the IEEE, 77, 257-286.}

\author{Jared O'Connell jared@well.ox.ac.uk}
\seealso{\code{hmmfit,hmmspec}}
\examples{
##See examples in 'hmmfit'
}