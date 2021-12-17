


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

set.seed(12434212)
#a = k and b = 1/lambda of the exponential
#Survival Model for exponential distributed data
n <- 10000
#mean survival = 1/lambda
beta <- c(5)
x <- cbind(rnorm(n,0,1))
lambda <- exp(beta*x)

time <- rexp(n,lambda)
data <- data.frame(time = time, x1 = x)

#now I constructed a data set and now estimate beta1-3
loglik_exp <- function(beta){
  #which is the best beta that leads to the best lambdas.
  -sum(dexp(data$time,exp(beta*data$x1), log = T))
}
optim(c(6),loglik_exp,method = "Brent",lower = -10,upper = 100)
-loglik_exp(beta)

#funktioniert im einidimensionalen
#daher jetzt zweidimensional:



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

set.seed(12434212)
#a = k and b = 1/lambda of the exponential
#Survival Model for exponential distributed data
n <- 10000
#mean survival = 1/lambda
beta <- c(-1,5)
x <- cbind(rnorm(n,0,1))
lambda <- exp(beta[1] + beta[2]*x)

time <- rexp(n,lambda)
data <- data.frame(time = time, x1 = x)

#now I constructed a data set and now estimate beta1-3
loglik_exp <- function(beta){
  #which is the best beta that leads to the best lambdas.
  -sum(dexp(data$time,exp(beta[1] + beta[2]*data$x1), log = T))
}
optim(c(0,0),loglik_exp)
-loglik_exp(beta)
#funktioniert auch im zweidimensionalen



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

set.seed(12434212)
#a = k and b = 1/lambda of the exponential
#Survival Model for exponential distributed data
n <- 100000
#mean survival = 1/lambda
beta <- c(-1,5,30)
x <- cbind(rnorm(n,0,1), rnorm(n,0,1))
lambda <- exp(beta[1] + beta[2]*x[,1] + beta[3]*x[,2])

#hier kommt die randomnes rein und deswegen muss die Lösung nicht gleich
#dkern Parametern sein
time <- rexp(n,lambda)
data <- data.frame(time = time, x1 = x[,1], x2 = x[,2])

#now I constructed a data set and now estimate beta1-3
loglik_exp <- function(beta){
  #which is the best beta that leads to the best lambdas.
  -sum(dexp(data$time,exp(beta[1] + beta[2]*data$x1+ beta[3]*data$x2), log = T))
}
#hier kommt es schon sehr auf die initialisierung an, damit das Maximum von der 
#Optim Funktion wirklich gefunden wird
optim(c(0,0,0),loglik_exp)
#wenn man mit den wahren Werten initiiert dann klappt es
optim(beta,loglik_exp)
#das Maximum ist gesucht
-loglik_exp(beta)
-loglik_exp(c(-11,2,16))
#funktioniert auch im zweidimensionalen



# 
# 
# #jetzt dreidimensional
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# 
# set.seed(12434212)
# #a = k and b = 1/lambda of the exponential
# #Survival Model for exponential distributed data
# beta <- c(0.1,0.5,10)
# n <- 10000
# x <- cbind(rnorm(n,0,1),rnorm(n,0,1))
# #mean survival = 1/lambda
# lambda <- exp(beta[1] + beta[2]*x[,1] + beta[3]*x[,2])
# 
# time <- rexp(n,lambda)
# data <- data.frame(time = time, x1 = x[,1], x2 = x[,2])
# 
# #now I constructed a data set and now estimate beta1-3
# loglik_exp <- function(beta){
#   -sum(dexp(data$time,exp(beta[1] + beta[2]*data$x1 + beta[3]*data$x2), log = T))
# }
# optim(c(6,6,6),loglik_exp)
# -loglik_exp(beta)
# #funktioniert nur so einigermaßen obwohl es keinen Zufall gibt, das Minimum
# #liegt nicht an der angepeilten Stelle