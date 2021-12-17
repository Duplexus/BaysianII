setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(coda)
library(ggplot2)
library(R2OpenBUGS)
####DATA SPECIFICATION####
data <- rweibull(100,2,1)
x <- sample(c(1:100)*0.01,100,T)

#Modell gegenüber wikipedia in openBugs: k=k lambda =1/lambda

# MODEL SPECIFICATION 
model.function <- function(){
  for (i in 1:N){
    y[i] ~ dweib(k, invlambda)
    log(taub[i]) <- beta0 * x[i]
    #log(invlambda[i]) <- beta0 * x[i]
  }
  #priors
  k ~ dunif(0.1,100)
  
  beta0 ~ dunif(-2,2)
  invlambda ~ dunif(-100,100)
}
write.model(model.function, "txt_files\\Weibull_test1.txt")
#all the stuff that is distributed
model.inits <- function(){list(k=2, beta0=1, invlambda = 2)}
#which data does the model need
model.data <- list( y = data, N = length(data), x = x)
# these are the parameters to save
parameters = c("k", "invlambda", "beta0", "taub")
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
