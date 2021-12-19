#Bayes1_Weibull_random.r
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(coda)
library(R2OpenBUGS)

Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")

model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group, id = Grub$id, Nsubj = length(unique(Grub$id)))

# MODEL SPECIFICATION 
model.function <- function(){
  for (i in 1:N){
    y[i] ~ dweib(k, invlambda[i])
    log(invlambda[i]) <- beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]]
  }
  #priors
  k ~ dunif(0.1,100)
  beta0 ~ dunif(-5,5)
  beta1 ~ dunif(-5,5)
  beta2 ~ dunif(-5,5)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,0.1)
  }
}
write.model(model.function, "Scripts\\Bayes1_Weibull_Random.txt")
model.inits <- function(){list(k=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20)) )}
parameters = c("k", "beta2", "beta0", "beta1", "b0")


model.out <- bugs(model.data, model.inits, 
                  model.file = "Bayes1_Weibull_Random.txt",
                  parameters=parameters,
                  n.chains = 2, n.iter = 5000,n.thin = 100,  n.burnin = 2000, debug = T,
                  codaPkg=T,
                  working.directory = ".\\Scripts")

Weibull_bayes <- read.bugs(model.out)

summary_weibull <- summary(Weibull_bayes)
b_params <- summary_weibull$statistics[grepl("b0\\[", dimnames(summary_weibull$statistics)[[1]]),]
mean(b_params[1:10,1])
attributes(summary_weibull$statistics)

