source("helpfunctions.r")
library("runjags")
library("coda")
library("dplyr")
Grub <- read.csv("../data/Grubs_Easy_normalized_size.csv")
Grub$number <- 1:nrow(Grub)
Grub <- Grub %>% arrange(upperlim)
length_Upper <- length(sort(Grub$upperlim))
#10.
lenngth_NA_Upper <- nrow(Grub) - length_Upper
NAs <- is.na(Grub$upperlim)
#just for numerical reasons, second has to be bigger than first
#and just to lazy to separate therefore makes no difference
Grub$upperlim[NAs] <- 12.000001
Grub$state <- c(rep(1,length_Upper),rep(2,lenngth_NA_Upper))
Grub$value2 <- as.numeric(NA)

model.inits <- list(list(thau2=2,sigma2_b0 =2 , beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20))),
                    list(thau2=20,sigma2_b0 =20, beta0=10, beta1 = 10,beta2 = 10 , b0 =  rnorm(20,0,30) ),
                    list(thau2=15,sigma2_b0 =15, beta0=20, beta1 = -10,beta2 = -15 , b0 =  runif(20,0,10))
)

# model.inits <-list(list(tau=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20))),
#                     list(tau=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20))))

parameters = c("sigma", "beta2", "beta0", "beta1", "b0", "Deviance","tau","y","thau2")

model.data <- list(y = Grub$value2, z = Grub$state, N1 = length_Upper,N2 = lenngth_NA_Upper,
                   x1 = Grub$grubsize, x2 = Grub$group, id = Grub$id,
                   lims = cbind(Grub$lowerlim,Grub$upperlim),
                   Nsubj = length(unique(Grub$id) ))

"
for (i in 1:N1){
  z[i] ~ dinterval(y[i], lims[i,])
  y[i] ~ dlnorm(mu[i], tau)
  mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
}
for (i in (N1+1):(N1+N2)){
  z[i] ~ dinterval(y[i], lims[i,])
  y[i] ~ dlnorm(mu[i], tau)
  mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
}
#priors
tau ~ dgamma(0.001, 0.001)
sigma <- sqrt(1/tau)

tau_b0 <- 1/sigma_b0
sigma_b0 ~ dunif(0,100)
beta0 ~ dnorm(0,10e-6)
beta1 ~ dnorm(0,10e-6)
beta2 ~ dnorm(0,10e-6)
for ( i in 1:Nsubj){
  b0[i] ~  dnorm(0,tau_b0)
}
}"
#### 2 Bigger beta priors ####
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  Deviance <- sum(D[])
  #priors
  thau2 ~ dunif(0,1000)
  tau <- 1/thau2
  sigma <- sqrt(1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  beta0 ~ dnorm(0,1e-3)
  beta1 ~ dnorm(0,1e-3)
  beta2 ~ dnorm(0,1e-3)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_cens_sens2 <- run.jags(model = model.function,
                               monitor = parameters, data = model.data,
                               inits = model.inits, burnin = 2000,
                               sample = 5000, thin = 1, n.chains = 3)
lognorm_rand_cens_sens_mcmc2 <- as.mcmc.list(lognorm_rand_cens_sens2)

##### 3 smaller beta priors ####
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  Deviance <- sum(D[])
  #priors
  thau2 ~ dunif(0,1000)
  tau <- 1/thau2
  sigma <- sqrt(1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  beta0 ~ dnorm(0,1e-9)
  beta1 ~ dnorm(0,1e-9)
  beta2 ~ dnorm(0,1e-9)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_cens_sens3 <- run.jags(model = model.function,
                               monitor = parameters, data = model.data,
                               inits = model.inits, burnin = 2000,
                               sample = 5000, thin = 1, n.chains = 3)
lognorm_rand_cens_sens_mcmc3 <- as.mcmc.list(lognorm_rand_cens_sens3)




####1 Vary the Distribution of the sigma and random effect sigma####
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  Deviance <- sum(D[])
  #priors
  thau2 ~ dunif(0,1000)
  tau <- 1/thau2
  sigma <- sqrt(1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  beta0 ~ dnorm(0,1e-06)
  beta1 ~ dnorm(0,1e-06)
  beta2 ~ dnorm(0,1e-06)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_cens_sens1 <- run.jags(model = model.function,
                               monitor = parameters, data = model.data,
                               inits = model.inits, burnin = 2000,
                               sample = 5000, thin = 1, n.chains = 3)
lognorm_rand_cens_sens_mcmc1 <- as.mcmc.list(lognorm_rand_cens_sens1)


####4 Vary the Distribution of the sigma and random effect sigma####
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  Deviance <- sum(D[])
  #priors
  thau2 ~ dgamma(0.2,0.2)
  tau <- thau2
  sigma <- sqrt(1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  beta0 ~ dnorm(0,1e-06)
  beta1 ~ dnorm(0,1e-06)
  beta2 ~ dnorm(0,1e-06)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_cens_sens4 <- run.jags(model = model.function,
                               monitor = parameters, data = model.data,
                               inits = model.inits, burnin = 2000,
                               sample = 5000, thin = 1, n.chains = 3)
lognorm_rand_cens_sens_mcmc4 <- as.mcmc.list(lognorm_rand_cens_sens4)


####5 Vary the Distribution of the random effect sigma####
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  Deviance <- sum(D[])
  #priors
  thau2 ~ dgamma(0.0001,0.0001)
  tau <- thau2
  sigma <- sqrt(1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dgamma(0.01,0.01)
  beta0 ~ dnorm(0,1e-06)
  beta1 ~ dnorm(0,1e-06)
  beta2 ~ dnorm(0,1e-06)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_cens_sens5 <- run.jags(model = model.function,
                               monitor = parameters, data = model.data,
                               inits = model.inits, burnin = 2000,
                               sample = 5000, thin = 1, n.chains = 3)


####6 Va everything####
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  Deviance <- sum(D[])
  #priors
  thau2 ~ dunif(0,1000)
  tau <- 1/thau2
  sigma <- sqrt(1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dgamma(0.01,0.01)
  beta0 ~ dnorm(0,1e-9)
  beta1 ~ dnorm(0,1e-9)
  beta2 ~ dnorm(0,1e-9)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_cens_sens6 <- run.jags(model = model.function,
                               monitor = parameters, data = model.data,
                               inits = model.inits, burnin = 2000,
                               sample = 5000, thin = 1, n.chains = 3)


#### 7 No Grubs ####
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta2 *x2[i]+ b0[id[i]]
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  Deviance <- sum(D[])
  #priors
  thau2 ~ dgamma(0.001,0.001)
  tau <- thau2
  sigma <- sqrt(1/tau)
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  beta0 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

runjags.options(method = "rjparallel")

lognorm_rand_cens_sens7 <- run.jags(model = model.function,
                               monitor = parameters, data = model.data,
                               inits = model.inits, burnin = 2000,
                               sample = 5000, thin = 1, n.chains = 3)
lognorm_rand_cens_sens_mcmc7 <- as.mcmc.list(lognorm_rand_cens_sens7)




