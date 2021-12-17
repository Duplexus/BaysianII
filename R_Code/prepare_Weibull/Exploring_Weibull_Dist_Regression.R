setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

set.seed(12434212)
#a = k and b = 1/lambda of the exponential
#Survival Model for Weibull distributed data
n <- 10000
#mean survival = 1/lambda
beta <- c(5)
x <- cbind(rnorm(n,0,1))
invlambda <- exp(beta*x)
k <- 2


time <- rweibull(n,k,invlambda)
data <- data.frame(time = time, x1 = x)

#now I constructed a data set and now estimate beta1-3
loglik_exp <- function(beta){
  #which is the best beta that leads to the best lambdas.
  -sum(dweibull(data$time,k,exp(beta*data$x1), log = T))
}
optim(c(6),loglik_exp,method = "Brent",lower = -10,upper = 100)
-loglik_exp(beta)
#funktioniert im einidimensionalen wenn ich invlambda modelle


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

set.seed(12434212)
#a = k and b = 1/lambda of the exponential
#Survival Model for exponential distributed data
n <- 10000
#mean survival = 1/lambda
beta <- c(5)
x <- cbind(rnorm(n,0,0.01))
invlambda <- 2
k <- exp(beta*x)


time <- rweibull(n,k,invlambda)
data <- data.frame(time = time, x1 = x)

#now I constructed a data set and now estimate beta1-3
loglik_exp <- function(beta){
  #which is the best beta that leads to the best lambdas.
  -sum(dweibull(data$time,exp(beta*data$x1),2, log = T))
}
optim(c(6),loglik_exp,method = "Brent",lower = -10,upper = 100)
loglik_exp(beta)
#funktioniert schlechter wenn man das k modelliert



