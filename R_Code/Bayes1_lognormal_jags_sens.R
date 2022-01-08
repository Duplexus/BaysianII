setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("helpfunctions.r")
library("runjags")
library("coda")
library("rjags")
Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")
#attributes(Grub$grubsize) <- NULL


#DEFINE INTITIAL VALUES
model.inits <- list(list(tau2=2, beta0=1, beta1 = 1,beta2 = 1 ),
                    list(tau2=20, beta0=10, beta1 = 10,beta2 = 10 ),
                    list(tau2=15, beta0=20, beta1 = -10,beta2 = -15 )
                    )

#Monitored Variables
parameters <-c("beta0", "beta1", "beta2", "sigma","tau")



#what happens if big grub size is forbidden



#Initial Values
model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group)

# Specification data model
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  #priors
  tau2 ~ dunif(0,1000)
  tau <- 1/tau2
  sigma <- sqrt(1/tau)
  beta0 ~ dnorm(0,10e-9)
  beta1 ~ dnorm(0,10e-9)
  beta2 ~ dnorm(0,10e-9)
}"

runjags.options(method = "rjparallel")
lognorm_sens <- run.jags(model = model.function,
                     monitor = parameters, data = model.data,
                     inits = model.inits, burnin = 2000,
                     sample = 5000, thin = 1, n.chains = 3)
lognorm_mcmc_sens <- as.mcmc.list(lognorm_sens)



model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  #priors
  tau2 ~ dunif(0,1000)
  tau <- 1/tau2
  sigma <- sqrt(1/tau)
  beta0 ~ dnorm(0,1)
  beta1 ~ dnorm(0,1)
  beta2 ~ dnorm(0,1)
}"

runjags.options(method = "rjparallel")
lognorm_sens2 <- run.jags(model = model.function,
                         monitor = parameters, data = model.data,
                         inits = model.inits, burnin = 2000,
                         sample = 5000, thin = 1, n.chains = 3)
