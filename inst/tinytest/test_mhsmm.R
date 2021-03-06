library(tinytest)
library(mhsmm)

#HMM
set.seed(1982)
J<-3
initial <- rep(1/J,J)
P <- matrix(c(.8,.5,.1,0.05,.2,.5,.15,.3,.4),nrow=J)
b <- list(mu=c(-10,0,10),sigma=c(2,1,.5))
model <- hmmspec(init=initial, trans=P, parms.emission=b,dens.emission=dnorm.hsmm)
model

train <- simulate(model, nsim=30, seed=1234, rand.emis=rnorm.hsmm)

init0 <- rep(1/J,J)
P0 <- matrix(1/J,nrow=J,ncol=J)
b0 <- list(mu=c(-5,0,5),sigma=c(1,1,1))
startval <- hmmspec(init=init0, trans=P0,parms.emission=b0,dens.emission=dnorm.hsmm)
h1 = hmmfit(train,startval,mstep=mstep.norm)

expect_equal(sum(train$s!=predict(h1,train)$s),0)

                                        #HSMM
set.seed(1982)
J <- 3
init <- c(0,0,1)
P <- matrix(c(0,1,0,0,0,1,1,0,0),nrow=J,byrow=TRUE)
B <- list(mu=c(-10,0,10),sigma=c(2,1,.5))
d <- list(lambda=c(10,30,60),shift=c(10,100,30),type='poisson')
model <- hsmmspec(init,P,parms.emission=B,sojourn=d,dens.emission=dnorm.hsmm)
train <- simulate(model,r=rnorm.hsmm,nsim=100,seed=123456)

start.poisson <- hsmmspec(init=rep(1/J,J)
                          ,transition=P
                          ,parms.emission=B
                          ,sojourn=list(lambda=c(9,25,40),shift=c(5,95,45),type='poisson')
                          ,dens.emission=dnorm.hsmm)

h.poisson <- hsmmfit(train,start.poisson,mstep=mstep.norm)

predicted <- predict(h.poisson,train)
expect_equal(sum(predicted$s!=train$s),0) ##there should be no errors

                                        #non-parametric
J <- 3
init <- c(0,0,1)
P <- matrix(c(0,1,0,0,0,1,1,0,0),nrow=J,byrow=TRUE)
B <- list(mu=c(-10,0,10),sigma=c(2,1,.5))

d=list(shape=c(100,200,400), scale = c(0.5,.75,0.25), type = "gamma")
model <- hsmmspec(init,P,parms.emission=B,sojourn=d,dens.emission=dnorm.hsmm)
train <- simulate.hsmmspec(model,r=rnorm.hsmm,nsim=100,seed=123456)

M <- 1:1000
d <- list(d=sapply(1:3,function(x) dgamma(M,model$sojourn$shape[x],1/model$sojourn$scale[x]))
         ,type='nonparametric')
model.start <- hsmmspec(init,P,parms.emission=B,sojourn=d,dens.emission=dnorm.hsmm)

h.np <- hsmmfit(train,model.start,mstep=mstep.norm)

predicted <- predict(h.np,train)
table(predicted$s,train$s)
expect_equal(sum(predicted$s!=train$s),0) ##there should be no errors

                                        # gamma sojourn

d=list(shape=c(100,200,400), scale = c(0.5,.75,0.25), type = "gamma")
model <- hsmmspec(init,P,parms.emission=B,sojourn=d,dens.emission=dnorm.hsmm)

train <- simulate(model,r=rnorm.hsmm,nsim=100,seed=123456)

h.gamma <- hsmmfit(train,model,mstep=mstep.norm)
predicted <- predict(h.gamma,train)
table(predicted$s,train$s)
expect_equal(sum(predicted$s!=train$s),0) ##there should be no errors


plot.sojourn(h.gamma$model)
plot.sojourn(model)
                                        #check it stops recurrent states
P <- matrix(c(0,1,0,0,0,1,0,0,1),nrow=J,byrow=TRUE)
expect_error(hsmmspec(init=rep(1/J,J)
                          ,transition=P
                          ,parms.emission=B
                          ,sojourn=list(lambda=c(9,25,40),shift=c(5,95,45),type='poisson')
                          ,dens.emission=dnorm.hsmm))

