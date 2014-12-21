\name{addStates}
\alias{addStates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Adds a bar representing state sequence. }
\description{
  Add a colour coded horizontal bar representing the state sequence to a plot of (presumably time-series) data.  
}
\usage{
addStates(states, x=NULL,ybot = axTicks(2)[1],
          ytop=ybot + (axTicks(2)[2] - axTicks(2)[1])/5,dy  = ytop - ybot,
          greyscale = FALSE,leg = NA, J = length(unique(states)), time.scale = 1, 
          shiftx = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{ 
  \item{states}{A vector of integers representing the states traversed}  
  \item{x}{The time values where the states are observed ((1:length(states)-shiftx)/time.scale if NULL)}
  \item{ybot}{Vertical bottom limit of the bar.}  
  \item{ytop}{Vertical top limit of the bar.}  
  \item{dy}{Height of the bar.}  
  \item{greyscale}{If TRUE produces a bar in greyscale.}  
  \item{leg}{Array of state names, if present, produces a legend.}  
  \item{J}{Number of states}  
  \item{time.scale}{Resolution of the timescale}  
  \item{shiftx}{Shift the bar forward or backwards horizontal by shiftx distance.}  
}
\seealso{ \code{addStates}}
\examples{
  plot(rnorm(100),type='l')
  addStates(rep(c(1,2),each=50))  

  plot(seq(0.01,1,.01),rnorm(100),type='l')
  addStates(rep(c(1,2),each=50),seq(0.01,1,.01))  
}
\author{Soren Hojsgaard sorenh@math.aau.dk}
