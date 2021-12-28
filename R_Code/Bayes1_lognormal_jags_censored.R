setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("helpfunctions.r")
library("runjags")
library("coda")
library("dplyr")
Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")
Grub <- Grub %>% arrange(upperlim)
length_Upper <- length(sort(Grub$upperlim))
#10
lenngth_NA_Upper <- nrow(Grub) - length_Upper
NAs <- is.na(Grub$upperlim)
#just for numerical reasons, upper has to be bigger than lower
Grub$upperlim[NAs] <- 12.000001
# PPO 
#da alle Intervallcensored sind gerade muss ich immer eine 1 schicken
#https://stats.stackexchange.com/questions/13847/how-does-dinterval-for-interval-censored-data-work-in-jags
Grub$state <- c(rep(1,length_Upper),rep(2,lenngth_NA_Upper))
#also helpfull to understand later syntax: https://stats.stackexchange.com/questions/70858/right-censored-survival-fit-with-jags
Grub$value <- as.numeric(NA)
model.data <- list(y = Grub$value, z = Grub$state, N1 = length_Upper,N2 = lenngth_NA_Upper, x1 = Grub$grubsize,
                    x2 = Grub$group, lims = cbind(Grub$lowerlim,Grub$upperlim))


model.inits <- list(list(sigma=2, beta0=1, beta1 = 1,beta2 = 1 ),list(sigma=2, beta0=1, beta1 = 1,beta2 = 1 ))
parameters <-c("beta0", "beta1", "beta2", "sigma")
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dlnorm(mu[i], sigma)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], sigma)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  #priors
  sigma ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
}"
runjags.options(method = "rjparallel")
ModelLogN <- run.jags(model = model.function,
                      monitor = parameters, data = model.data,
                      inits = model.inits, burnin = 1000, sample = 10000, thin = 1, n.chains = 2
                      , progress.bar = "text")


#Model Diagnostics and so on 
#Starts with ModelLogN and plots it and calculates BIC and the 
#Statistics of Rubin and so on
#some Results
plot(ModelLogN)
print(ModelLogN)
#The DIC Value for model comparison
dic_val <- extract.runjags(ModelLogN, "dic")
dic_val
extract.runjags(ModelLogN, "stochastic")

#coda integration so also coda stuff is available
#Model Diagnostik plots
mcmc <- as.mcmc.list(ModelLogN)
gelman.diag(mcmc, confidence = 0.95)
gelman.plot(mcmc, confidence = 0.95)
geweke.diag(mcmc)
geweke.plot(mcmc)
plot(mcmc)


