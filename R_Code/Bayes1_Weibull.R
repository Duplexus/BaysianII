#Weibull Modell 1 no random effects
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(coda)
library(R2OpenBUGS)
Grub <- read.csv("..\\data\\Grubs_Easy.csv")
model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group, )

# MODEL SPECIFICATION 
model.function <- function(){
  for (i in 1:N){
    y[i] ~ dweib(k, invlambda[i])
    log(invlambda[i]) <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  #priors
  k ~ dunif(0.1,100)
  beta0 ~ dunif(-5,5)
  beta1 ~ dunif(-5,5)
  beta2 ~ dunif(-5,5)
}
write.model(model.function, "Scripts\\Bayes1_Weibull.txt")
model.inits <- function(){list(k=2, beta0=1, beta1 = 1,beta2 = 1 )}
parameters = c("k", "beta2", "beta0", "beta1")


model.out <- bugs(model.data, model.inits, 
                  model.file = "Bayes1_Weibull.txt",
                  parameters=parameters,
                  n.chains = 3, n.iter = 5000,  n.burnin = 1000, debug = T,
                  codaPkg=T,
                  working.directory = ".\\Scripts")

Weibull_bayes <- read.bugs(model.out)

summary(Weibull_bayes)

