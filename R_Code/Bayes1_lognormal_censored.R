setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
Grub <- read.csv("..\\data\\Grubs_Nematodes.csv")
library(survival)
library(coda)
library(R2OpenBUGS)
NAs <- is.na(Grub$UPPERLIM)
Grub$EVENT <- ifelse(is.na(Grub$UPPERLIM),0,3)
Grub$UPPERLIM[NAs]<- 14
Grub$VALUE <- (Grub$UPPERLIM + Grub$LOWERLIM) /2
# Grub$VALUE  <- as.numeric(NA)

y <- unlist(as.vector(Grub$VALUE))
y <-Grub$VALUE
#y[1:1] <- NA
#y[81:140] <- NA

model.data <- list(lower = Grub$LOWERLIM, upper = Grub$UPPERLIM, N = length(Grub$UPPER), x1 = Grub$GRUBSIZE,
                    x2 = Grub$GROUP, y = y)



# MODEL SPECIFICATION 
model.function <- function(){
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], sigma)C(lower[i], upper[i])
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  #priors
  sigma ~ dunif(0.1,100)
  beta0 ~ dunif(-5,5)
  beta1 ~ dunif(-5,5)
  beta2 ~ dunif(-5,5)
}
#write.model(model.function, "Scripts\\Bayes1_lognormal_censored.txt")
model.inits <- function(){list(sigma=2, beta0=1, beta1 = 1,beta2 = 1)}
parameters = c("sigma", "beta2", "beta0", "beta1")


model.out <- bugs(model.data, model.inits, 
                  model.file = "Bayes1_lognormal_censored.txt",
                  parameters=parameters,
                  n.chains = 1, n.iter = 5000,  n.burnin = 2000, debug = F,
                  codaPkg=T,
                  working.directory = ".\\Scripts")
lognormal_bayes2 <- read.bugs(model.out)

summary(lognormal_bayes2)
summary(lognormal_bayes)






















# 
# lognormal_results <- survreg(Surv(time = LOWERLIM,time2 = UPPERLIM, event = EVENT) ~ GRUBSIZE +  GROUP, 
#                              Grub, dist = "lognormal")
# summary(lognormal_results)
# 
# Surv(Grub$LOWERLIM,Grub$UPPERLIM, event = Grub$EVENT)
# require(flexsurv)
# flexsurvreg(formula = Surv(time = LOWERLIM, time2 = UPPERLIM, event = EVENT) ~ GRUBSIZE + GROUP, 
#             data=Grub, dist = "lognormal") 









