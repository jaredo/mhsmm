\name{predict.hsmm}
\alias{predict.hsmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction for hsmms}
\description{
  Predicts the underlying state sequence for an observed sequence \code{newdata} given a \code{hsmm} model
}
\usage{
\method{predict}{hsmm}(object, newdata, method = "viterbi", ...)
}

\arguments{
  \item{object}{ An object of type \code{hsmm}}
  \item{newdata}{ A vector or dataframe of observations }
  \item{method}{Prediction method (see details)}
  \item{\dots}{further arguments passed to or from other methods.}

}
\details{
    If \code{method="viterbi"}, this technique applies the Viterbi algorithm for HSMMs, producing the most likely sequence of states given the observed data.  If \code{method="smoothed"}, then the individually most likely (or smoothed) state sequence is produced, along with a matrix with the respective probabilities for each state.
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

\seealso{ \code{\link{hsmmfit}},\code{\link{predict.hsmmspec}}}
\examples{
##See 'hsmmfit' for examples
}