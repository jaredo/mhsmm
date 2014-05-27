\name{gammafit}
\alias{gammafit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parameter estimation for the Gamma distribution}
\description{
Estimates parameters for the Gamma distribution using the Method of Maximum Likelihood, works with weighted data.
}
\usage{
gammafit(x, wt = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of observations}
  \item{wt}{Optional set of weights}
}
\value{
  \item{shape}{The shape parameter}
  \item{scale}{The scale parameter (equal to 1/rate)}
}
\references{Choi, S. and Wette, R. (1969), Maximum likelihood estimation of the parameters of the gamma distribution and their bias, Technometrics, 11, 683-96-690.}
\author{Jared O'Connell jaredoconnell@gmail.com}
\examples{
  gammafit(rgamma(1000,shape=10,scale=13))
}
