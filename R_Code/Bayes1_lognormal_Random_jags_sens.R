#Bayes1_lognormal_random.r
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("runjags")
library("coda")
library("rjags")

Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")

Grub$grubsize
# 
# model.inits <-list(list(tau=2, beta1 = 1, b0 = c(rep(1,times = 20))),
#                    list(tau=2, beta1 = 1, b0 = c(rep(1,times = 20))))


model.inits <- list(list(tau2=2, beta0=1, beta1 = 1,beta2 = 1 ),
                    list(tau2=20, beta0=10, beta1 = 10,beta2 = 10 ),
                    list(tau2=15, beta0=20, beta1 = -10,beta2 = -15 )
)
parameters = c("sigma2", "beta2", "beta1","beta0","sigma2_b0","tau")
#parameters = c("sigma", "beta1", "b0")

model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    id = Grub$id,x2 = Grub$group,  Nsubj = length(unique(Grub$id)))
#### 2 Bigger beta priors ####
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    
    
  }
  #priors
  tau2 ~ dunif(0,1000)
  tau <- 1/tau2
  sigma2 <- (1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_sens2 <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 2000,
                    sample = 5000, thin = 1, n.chains = 3)
lognorm_rand_sens_mcmc2 <- as.mcmc.list(lognorm_rand_sens2)
dic2 <- extract.runjags(lognorm_rand_sens2, "dic")

# summary(lognorm_rand_sens2)
##### 3 smaller beta priors ####
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    
    
  }
  #priors
  tau2 ~ dunif(0,1000)
  tau <- 1/tau2
  sigma2 <- (1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  beta0 ~ dnorm(0,10e-9)
  beta1 ~ dnorm(0,10e-9)
  beta2 ~ dnorm(0,10e-9)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_sens3 <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 2000,
                    sample = 5000, thin = 1, n.chains = 3)
lognorm_rand_sens_mcmc3 <- as.mcmc.list(lognorm_rand_sens3)

summary(lognorm_rand_sens3)
dic3 <- extract.runjags(lognorm_rand_sens3, "dic")



####1 Vary the Distribution of the sigma and random effect sigma####
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    
    
  }
  #priors
  tau2 ~ dunif(0,1000)
  tau <- 1/tau2
  sigma2 <- (1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_sens1 <- run.jags(model = model.function,
                               monitor = parameters, data = model.data,
                               inits = model.inits, burnin = 2000,
                               sample = 5000, thin = 1, n.chains = 3)
lognorm_rand_sens_mcmc1 <- as.mcmc.list(lognorm_rand_sens1)
dic1 <- extract.runjags(lognorm_rand_sens1, "dic")

summary(lognorm_rand_sens_mcmc1)
summary(lognorm_rand_sens1)

####4 Vary the Distribution of the sigma and random effect sigma####
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    
    
  }
  #priors
  tau2 ~ dgamma(0.2,0.2)
  tau <- tau2
  sigma2 <- (1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_sens4 <- run.jags(model = model.function,
                               monitor = parameters, data = model.data,
                               inits = model.inits, burnin = 2000,
                               sample = 5000, thin = 1, n.chains = 3)
lognorm_rand_sens_mcmc4 <- as.mcmc.list(lognorm_rand_sens4)
dic4 <- extract.runjags(lognorm_rand_sens4, "dic")

# 
# summary(lognorm_rand_sens_mcmc4)
# summary(lognorm_rand_sens4)

####5 Vary the Distribution of the random effect sigma####
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    
    
  }
  #priors
  tau2 ~ dgamma(0.0001,0.0001)
  tau <- tau2
  sigma2 <- (1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dgamma(0.01,0.01)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_sens5 <- run.jags(model = model.function,
                               monitor = parameters, data = model.data,
                               inits = model.inits, burnin = 2000,
                               sample = 5000, thin = 1, n.chains = 3)

# summary(lognorm_rand_sens5)
dic5 <- extract.runjags(lognorm_rand_sens5, "dic")

####6 Va everything####
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0+ beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    
    
  }
  #priors
  tau2 ~ dunif(0,1000)
  tau <- 1/tau2
  sigma2 <- (1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dgamma(0.01,0.01)
  beta0 ~ dnorm(0,10e-9)
  beta1 ~ dnorm(0,10e-9)
  beta2 ~ dnorm(0,10e-9)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_sens6 <- run.jags(model = model.function,
                               monitor = parameters, data = model.data,
                               inits = model.inits, burnin = 2000,
                               sample = 5000, thin = 1, n.chains = 3)

summary(lognorm_rand_sens6)
dic6 <- extract.runjags(lognorm_rand_sens6, "dic")

#### 7 No Grubs ####
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta2 *x2[i]+ b0[id[i]]
    
    
  }
  #priors
  tau2 ~ dgamma(0.001,0.001)
  tau <- tau2
  sigma2 <- (1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  beta0 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_sens7 <- run.jags(model = model.function,
                               monitor = parameters, data = model.data,
                               inits = model.inits, burnin = 2000,
                               sample = 5000, thin = 1, n.chains = 3)
lognorm_rand_sens_mcmc7 <- as.mcmc.list(lognorm_rand_sens7)

summary(lognorm_rand_sens7)
dic7 <- extract.runjags(lognorm_rand_sens7, "dic")
