\name{smooth.discrete}
\alias{smooth.discrete}
\alias{print.smoothDiscrete}
\alias{summary.smoothDiscrete}
\alias{predict.smoothDiscrete}
\alias{createTransition}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Smoothing a discrete time series.
}
\description{
  The smooth.discrete() function provides a simple smoothing of a
  time series of discrete values measured at equidistant times.
  Under the hood of  smooth.discrete() is a hidden Markov model.
}
\usage{
smooth.discrete(y, init = NULL, trans = NULL, parms.emission = 0.5, 
                method = "viterbi", details = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    A numeric vector
  }
  \item{init}{
    Initial distribution (by default derived from data; see the vignette
    for details)
}
  \item{trans}{
    Transition matrix (by default derived from data; see the vignette
    for details)
}
\item{parms.emission}{
  Matrix describing the conditional
  probabilities of the observed states given the latent states. (See the
  vignette for details).
}
  \item{method}{
    Either "viterbi" or "smoothed". The viterbi method gives the jointly
    most likely sequence; the smoothed method gives the sequence of
    individually most likely states.
}
  \item{details}{
    Controlling the amount of information printed.
}
  \item{\dots}{
    Further arguments passed on to the "hmmfit" function.
}
}
\details{
 The parameters are estimated using the Baum-Welch algorithm (a special
 case of the EM-algorithm).
}
\value{
  A list with the following components:
  \item{s}{The "smoothed" states}
  \item{model}{The underlying hmm (hidden Markov model) object}
  \item{data}{The data}
  \item{initial}{The initial parameters}
}
%\references{  
%}
\author{
  Søren Højsgaard <sorenh at agrsci.dk>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hmmspec}}, \code{\link{hmmfit}}
}
\examples{
## Please see the vignette
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

