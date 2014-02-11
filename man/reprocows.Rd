\name{reprocows}
\docType{data}
\alias{reprocows}
\title{Reproductive data from seven dairy cows}
\description{
  This data set contains hourly observations on progesterone and an activity index at hourly intervals since calving on seven dairy cows.  
}
\usage{reprocows}
\format{\code{reprocows} is a data frame containing 13040 rows. 
\code{id} is the cow ID, \code{progesterone} is a measurement of the
hormone in ng/L taken from a milk sample, activity is a relative measure
of activity calculated from a pedometer. 

There are a large number of missing values as progesterone is measured
only at milking time (and at a farm manager's discretion).  Missing
values in activity occur due to hardware problems can occur with
pedometers. 
}
\source{Danish Cattle Research Centre}
\references{
Peters, A. and Ball, P. (1995), "Reproduction in Cattle," 2nd ed.
}
\keyword{datasets}

\examples{
data(reprocows)
data(reproai)
data(reproppa)
tm = 1600

J <- 3
init <- c(1,0,0)
trans <- matrix(c(0,0,0,1,0,1,0,1,0),nrow=J)
emis <- list(mu=c(0,2.5,0),sigma=c(1,1,1))

N <- as.numeric(table(reprocows$id))
train <- list(x=reprocows$activity,N=N)
class(train) <- "hsmm.data"
tmp <- gammafit(reproppa * 24)
M <- max(N)

d <- cbind(dgamma(1:M,shape=tmp$shape,scale=tmp$scale),
 # ppa sojourn directly estimated from ppa data set
dunif(1:M,4,30),
 # oestrus between 4 and 30 hours
dunif(1:M,15*24,40*24))
 #not-oestrus between 15 and 40 days

startval <- hsmmspec(init,trans,parms.emission=emis,list(d=d,type='gamma'),
  dens.emission=dnorm.hsmm)
#not run (takes some time)
#h.activity <- hsmmfit(train,startval,mstep=mstep.norm,maxit=10,M=M,lock.transition=TRUE)
  
}
