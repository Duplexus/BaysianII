#Bayes1_lognormal_random.r
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library("runjags")
library("coda")
library("rjags")

Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")

model.inits <-list(list(sigma=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20))),
                   list(sigma=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20))))

parameters = c("sigma", "beta2", "beta0", "beta1", "b0")

model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group, id = Grub$id, Nsubj = length(unique(Grub$id)))



model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], sigma)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
  }
  #priors
  sigma ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,0.0001)
  }
}"


lognorm_rand <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 2000,
                    sample = 5000, thin = 1, n.chains = 2)


parameters = c("sigma", "beta2", "beta0", "beta1", "b0","predict","b0.rep","ppo")
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], sigma)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    predict[i]  ~ dlnorm(mu[i], sigma)
  }
  # forecast <- y[] - predict[]
  #priors
  sigma ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,0.0001)
    # Distribution of future b0_i
    b0.rep[i] ~ dnorm(0,0.0001)  
  }
  for (i in 1:N){
    ppo[i] <- dlnorm(y[i],mu[i],sigma)
  }
}"


lognorm_rand_rep <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 2000,
                    sample = 5000, thin = 1, n.chains = 2)
lognorm_rand_mcmc_rep <- as.mcmc.list(lognorm_rand_rep)

