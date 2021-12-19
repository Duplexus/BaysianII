#lognormal Modell 1 no random effects
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(coda)
library(R2OpenBUGS)
Grub <- read.csv("..\\data\\Grubs_Easy.csv")
Grub$grubsize <- scale(Grub$grubsize)
attributes(Grub$grubsize) <- NULL
# Grub$group <- scale(Grub$group)

model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group )

# MODEL SPECIFICATION 
model.function <- function(){
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], sigma)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  #priors
  sigma ~ dunif(0.1,100)
  beta0 ~ dunif(-5,5)
  beta1 ~ dunif(-5,5)
  beta2 ~ dunif(-5,5)
}
write.model(model.function, "Scripts\\Bayes1_lognormal.txt")
model.inits <- function(){list(sigma=2, beta0=1, beta1 = 1,beta2 = 1 )}
parameters = c("sigma", "beta2", "beta0", "beta1")


model.out <- bugs(model.data, model.inits, 
                  model.file = "Bayes1_lognormal.txt",
                  parameters=parameters,
                  n.chains = 1, n.iter = 5000,  n.burnin = 1000, debug = T,
                  codaPkg=T,
                  working.directory = ".\\Scripts")

lognormal_bayes <- read.bugs(model.out)

summary(lognormal_bayes)

#parametric
library(survival)
lognormal_results <- survreg(Surv(value) ~ grubsize +  group, Grub, dist = "lognormal")
# scale is sigma seems to be 1/the one from the baysian but beside this similar
#results between bayes and model are really similar
summary(lognormal_results)



# but influence is estimated similar
1/0.62
y <- rweibull(1000, shape=2, scale=5)
survreg(Surv(y)~1, dist="weibull")


