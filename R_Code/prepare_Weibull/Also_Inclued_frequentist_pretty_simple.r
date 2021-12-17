setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

hist(rweibull(10000,1,1),breaks = 100)
hist(exp(rnorm(10000)),breaks = 100)
x <- 1:1000 * 0.002

#a = k and b = 1/lambda of the exponential
plot(x,dweibull(x,1,4), type = "l", ylim = c(0,2))
lines(x,dexp(x,0.5)+0.1, col = "red")
#both lines are the same (+0.1) that one can see

set.seed(12434212)
#a = k and b = 1/lambda of the exponential
data <- rweibull(100,10,1)

optim()

loglik_weib <- function(x){
  -sum(dweibull(data,x[1],x[2], log = T))
}
optim(c(6,6),loglik_weib)

#so kann man parametrisch eine Weibull-verteilung schätzen
#Das ganze jetzt mit Bayes schätzen:

#############################################
library(coda)
library(ggplot2)
library(R2OpenBUGS)

# MODEL SPECIFICATION 
model.function <- function(){
    for (i in 1:N){
      x[i] ~ dweib(k, invlambda)
          }
    #priors
    k ~ dunif(0,100)
    invlambda ~ dunif(0,100)
}
write.model(model.function, "txt_files\\Weibull_test1.txt")
model.inits <- function(){list(k=2, invlambda=2)}
model.data <- list( x = data, N = length(data))
# these are the parameters to save
parameters = c("k", "invlambda")
getwd()
# specify model, data, number of parallel chains
model.out <- bugs(model.data, model.inits, 
                  model.file = "Weibull_test1.txt",
                  parameters=parameters,
                  n.chains = 1, n.iter = 5000,  n.burnin = 100, debug = T,
                  codaPkg=T,
                  working.directory = ".\\txt_files")

Weibull_bayes <- read.bugs(model.out)
#in openbugs lambda is the same as for the exponential
#that means 1-exp(lambda*x) and not as in R 1-exp(x/lambda) 
summary(Weibull_bayes)

