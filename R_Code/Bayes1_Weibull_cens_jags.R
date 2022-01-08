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
Grub$state <- c(rep(1,length_Upper),rep(2,lenngth_NA_Upper))
Grub$value2 <- as.numeric(NA)
model.data <- list(y = Grub$value2, z = Grub$state, N1 = length_Upper,N2 = lenngth_NA_Upper, x1 = Grub$grubsize,
                   x2 = Grub$group, lims = cbind(Grub$lowerlim,Grub$upperlim))


# MODEL SPECIFICATION 
#now the model specification is in line with the estimation in survival
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  #priors
  scale <- 1/k
  k ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
}"
model.inits <- list(list(k=3, beta0=1, beta1 = 1,beta2 = 1 ),
                    list(k=1, beta0=10, beta1 = 10,beta2 = 10 ),
                    list(k=0.5, beta0=20, beta1 = -10,beta2 = -15 )
)

#Monitored Variables
parameters <-c("beta0", "beta1", "beta2","scale")

runjags.options(method = "rjparallel")
weibull_cens <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 2000,
                    sample = 5000, thin = 1, n.chains = 3)

weibull_cens_mcmc <- as.mcmc.list(weibull_cens)


model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dweib(k, invlambda[i])
    y_rep[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dweib(k, invlambda[i])
    y_rep[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  #priors
  scale <- 1/k
  k ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  for (i in 1:(N1+N2)){
    ppo[i] <- dweib(y[i],k, invlambda[i])
    ppo_rep[i] <- dweib(y_rep[i],k, invlambda[i])
    #for DIC
    D[i] <- -2*log(dweib(y[i],k, invlambda[i]))
  }
  Deviance <- sum(D[])
}"

parameters <-c("ppo_rep","y_rep","beta0", "beta1", "beta2", "invlambda","scale","ppo","Deviance","y")
weibull_cens_rep <- run.jags(model = model.function,
                          monitor = parameters, data = model.data,
                          inits = model.inits, burnin = 2000,
                          sample = 5000, thin = 1, n.chains = 3)

weibull_cens_rep_mcmc <- as.mcmc.list(weibull_cens_rep)
