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
#just for numerical reasons, second has to be bigger than first
#and just to lazy to separate therefore makes no difference
Grub$upperlim[NAs] <- 12.000001
Grub$state <- c(rep(1,length_Upper),rep(2,lenngth_NA_Upper))
Grub$value2 <- as.numeric(NA)

model.inits <-list(list(sigma=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20))),
                    list(sigma=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20))))

parameters = c("sigma", "beta2", "beta0", "beta1", "b0")

model.data <- list(y = Grub$value2, z = Grub$state, N1 = length_Upper,N2 = lenngth_NA_Upper,
                   x1 = Grub$grubsize, x2 = Grub$group, id = Grub$id,
                   lims = cbind(Grub$lowerlim,Grub$upperlim), N = length(Grub$value),
                   Nsubj = length(unique(Grub$id) ))


model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], sigma)
    #mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    mu[i] <- beta1 *x1[i] + b0[id[i]]
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], sigma)
    #mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    mu[i] <- beta1 *x1[i] + b0[id[i]]
  }
  #priors
  sigma ~ dgamma(0.001, 0.001)
  # beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  # beta2 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,0.0001)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_cens <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 2000,
                    sample = 5000, thin = 1, n.chains = 2)
lognorm_rand_cens_mcmc <- as.mcmc.list(lognorm_rand_cens)

parameters = c("sigma", "beta2", "beta0", "beta1", "b0","predict","b0.rep","ppo")
model.function <- "model{
   for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], sigma)
    #mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    mu[i] <- beta1 *x1[i] + b0[id[i]]
    predict[i]  ~ dlnorm(mu[i], sigma)
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], sigma)
    #mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    mu[i] <- beta1 *x1[i] + b0[id[i]]
    predict[i]  ~ dlnorm(mu[i], sigma)
  }
  
  
  
  
  # forecast <- y[] - predict[]
  #priors
  sigma ~ dgamma(0.001, 0.001)
  # beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  # beta2 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,0.0001)
    # Distribution of future b0_i
    b0.rep[i] ~ dnorm(0,0.0001)  
  }
  for (i in 1:N){
    ppo[i] <- dlnorm(y[i],mu[i],sigma)
  }
}"


lognorm_rand_cens_rep <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 2000,
                    sample = 5000, thin = 1, n.chains = 2)
lognorm_rand_cens_mcmc_rep <- as.mcmc.list(lognorm_rand_cens_rep)

