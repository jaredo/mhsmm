\name{hsmmfit}
\alias{hsmmfit}
\title{fit a hidden semi-Markov model}
\description{
  Estimates parameters of a HSMM using the EM algorithm.
}
\usage{
hsmmfit(x,model,mstep=NULL,M=NA,maxit=100,
        lock.transition=FALSE,lock.d=FALSE,graphical=FALSE)
}
\arguments{
  \item{x}{A hsmm.data object (see Details)}
  \item{model}{ Starting parameters for the model (see \code{hsmmspec}) }
  \item{mstep}{ Re-estimates the parameters of density function on each iteration }
  \item{maxit}{ Maximum number of iterations }
  \item{M}{Maximum number of time spent in a state (truncates the waiting distribution)} 
  \item{lock.transition}{ If TRUE will not re-estimate the transition matrix }
  \item{lock.d}{ If TRUE will not re-estimate the sojourn time density }
  \item{graphical}{If TRUE will plot the sojourn densities on each iteration}  
}
\value{
  \item{start}{A vector of the starting probabilities for each state}
  \item{a}{The transition matrix of the embedded Markov chain}
  \item{emission}{A list of the parameters of the emission distribution}
  \item{waiting}{A list of the parameters of the waiting distribution}	
}
\references{Jared O'Connell, Soren Hojsgaard (2011). Hidden Semi Markov Models for
Multiple Observation Sequences: The mhsmm Package for R., 
Journal of Statistical Software, 39(4), 1-22., URL
http://www.jstatsoft.org/v39/i04/.

Guedon, Y. (2003),  Estimating hidden semi-Markov chains from discrete sequences, 
	Journal of Computational and Graphical Statistics, Volume 12, Number 3, page 604-639 - 2003 }
\author{ Jared O'Connell jaredoconnell@gmail.com }
\seealso{\code{\link{hsmmspec}},\code{\link{simulate.hsmmspec}},\code{\link{predict.hsmm}}
}
\examples{
J <- 3
init <- c(0,0,1)
P <- matrix(c(0,.1,.4,.5,0,.6,.5,.9,0),nrow=J)
B <- list(mu=c(10,15,20),sigma=c(2,1,1.5))
d <- list(lambda=c(10,30,60),shift=c(10,100,30),type='poisson')
model <- hsmmspec(init,P,parms.emission=B,sojourn=d,dens.emission=dnorm.hsmm)
train <- simulate(model,r=rnorm.hsmm,nsim=100,seed=123456)
plot(train,xlim=c(0,400))
start.poisson <- hsmmspec(init=rep(1/J,J),
  transition=matrix(c(0,.5,.5,.5,0,.5,.5,.5,0),nrow=J),
  parms.emission=list(mu=c(4,12,23),
		sigma=c(1,1,1)),
  sojourn=list(lambda=c(9,25,40),shift=c(5,95,45),type='poisson'),
 dens.emission=dnorm.hsmm)

M=500
# not run (takes some time)
# h.poisson <- hsmmfit(train,start.poisson,mstep=mstep.norm,M=M)
# plot(h.poisson$loglik,type='b',ylab='Log-likelihood',xlab='Iteration') ##has it converged?
# summary(h.poisson)
# predicted <- predict(h.poisson,train)  
# table(train$s,predicted$s) ##classification matrix
# mean(predicted$s!=train$s) ##misclassification rate

d <- cbind(dunif(1:M,0,50),dunif(1:M,100,175),dunif(1:M,50,130))
start.np <- hsmmspec(init=rep(1/J,J),
  transition=matrix(c(0,.5,.5,.5,0,.5,.5,.5,0),nrow=J),
  parms.emission=list(mu=c(4,12,23),
  sigma=c(1,1,1)),
  sojourn=list(d=d,type='nonparametric'),
  dens.emission=dnorm.hsmm)
# not run (takes some time)
# h.np <- hsmmfit(train,start.np,mstep=mstep.norm,M=M,graphical=TRUE)
# mean(predicted$s!=train$s) ##misclassification rate


#J <- 2
#init <- c(1, 0)
#P <- matrix(c(0, 1, 1, 0), nrow = J)
#B <- list(mu = list(c(2, 3), c(3, 4)), sigma = list(matrix(c(4, 2, 2, 3), ncol = 2), diag(2)))
#d <- list(shape = c(10, 25), scale = c(2, 2), type = "gamma")
#model <- hsmmspec(init, P, parms.emis = B, sojourn = d, dens.emis = dmvnorm.hsmm)
#train <- simulate(model, c(1000,100), seed = 123, rand.emis = rmvnorm.hsmm)

#yhat <- predict(model, train)
#mean(predict(model,train)$s==train$s)


}

