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
# PPO 
#da alle Intervallcensored sind gerade muss ich immer eine 1 schicken
#https://stats.stackexchange.com/questions/13847/how-does-dinterval-for-interval-censored-data-work-in-jags
Grub$state <- c(rep(1,length_Upper),rep(2,lenngth_NA_Upper))
#also helpfull to understand later syntax: https://stats.stackexchange.com/questions/70858/right-censored-survival-fit-with-jags
Grub$value2 <- as.numeric(NA)
model.data <- list(y = Grub$value2, z = Grub$state, N1 = length_Upper,N2 = lenngth_NA_Upper, x1 = Grub$grubsize,
                    x2 = Grub$group, lims = cbind(Grub$lowerlim,Grub$upperlim))


model.inits <- list(list(tau=2, beta0=1, beta1 = 1,beta2 = 1 ),
                    list(tau=20, beta0=10, beta1 = 10,beta2 = 10 ),
                    list(tau=40, beta0=-100, beta1 = 100,beta2 = -30 )
)
parameters <-c("beta0", "beta1", "beta2", "tau","sigma2")
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
    #for DIC
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  #priors
  tau ~ dgamma(0.001, 0.001)
  sigma2 <- (1/tau)
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
}"
runjags.options(method = "rjparallel")
lognorm_cens <- run.jags(model = model.function,
                      monitor = parameters, data = model.data,
                      inits = model.inits, burnin = 1000, sample = 10000, thin = 1, n.chains = 3
                      , progress.bar = "text")
lognorm_cens_mcmc <- as.mcmc.list(lognorm_cens)

# #Model Diagnostics and so on 
# #Starts with ModelLogN and plots it and calculates BIC and the 
# #Statistics of Rubin and so on
# #some Results
# plot(ModelLogN)
# print(ModelLogN)
# #The DIC Value for model comparison
# dic_val <- extract.runjags(ModelLogN, "dic")
# dic_val
# extract.runjags(ModelLogN, "stochastic")
# 
# #coda integration so also coda stuff is available
# #Model Diagnostik plots
# mcmc <- as.mcmc.list(ModelLogN)
# gelman.diag(mcmc, confidence = 0.95)
# gelman.plot(mcmc, confidence = 0.95)
# geweke.diag(mcmc)
# geweke.plot(mcmc)
# plot(mcmc)

model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dlnorm(mu[i], tau)
    y_rep[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
    #for DIC
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], tau)
    y_rep[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
    #for DIC
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  Deviance <- sum(D[])
  #priors
  tau ~ dgamma(0.001, 0.001)
  sigma <- sqrt(1/tau)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
    #ppc look on residuals
  for (i in 1:(N1+N2)){
    ppo[i] <- dlnorm(y[i],mu[i],tau)
    ppo_rep[i] <- dlnorm(y_rep[i],mu[i],tau)
    k[i] <- log(y[i])
    res[i] <- k[i] - mu[i]
  }
}"
runjags.options(method = "rjparallel")
parameters <-c("y_rep","beta0", "beta1", "beta2", "tau","ppo","res","mu",
               "Deviance","y","sigma")
lognorm_cens_rep <- run.jags(model = model.function,
                         monitor = parameters, data = model.data,
                         inits = model.inits, burnin = 2000, sample = 5000, thin = 1, n.chains = 3
                         , progress.bar = "text")
lognorm_cens_rep_mcmc <- as.mcmc.list(lognorm_cens_rep)


