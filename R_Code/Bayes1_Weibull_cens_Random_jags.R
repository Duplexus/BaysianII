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

parameters = c("sigma", "beta2", "beta0", "beta1", "b0")

model.data <- list(y = Grub$value2, z = Grub$state, N1 = length_Upper,N2 = lenngth_NA_Upper,
                   x1 = Grub$grubsize, x2 = Grub$group, id = Grub$id,
                   lims = cbind(Grub$lowerlim,Grub$upperlim),
                   Nsubj = length(unique(Grub$id)))
# MODEL SPECIFICATION 
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dweib(k, invlambda[i])
    #log(invlambda[i]) <- beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]]
    log(invlambda[i]) <- beta1 *x1[i] + b0[id[i]]
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dweib(k, invlambda[i])
    #log(invlambda[i]) <- beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]]
    log(invlambda[i]) <- beta1 *x1[i] + b0[id[i]]
  }
  #priors
  k ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  sigma ~ dunif(0,100)
  beta2 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(mean(beta0 + beta2 *x2[i:(i+6)]) ,sigma)
    b0_rep[i] ~ dnorm(mean(beta0 + beta2 *x2[i:(i+6)]) ,sigma)
  }
}"

#Monitored Variables
parameters <- c("k", "sigma", "beta2", "beta0", "beta1","b0", "b0_rep")
#parameters <- c("k", "beta1", "b0", "b0_rep", "sigma")

model.inits <- list(list(k=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20)),b0_rep = c(rep(1,times = 20))),
                    list(k=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20)),b0_rep = c(rep(1,times = 20))))
# model.inits <- list(list(k=2,beta1 = 1, b0 = c(rep(1,times = 20))),
#                     list(k=2,beta1 = 1, b0 = c(rep(1,times = 20))))


runjags.options(method = "rjparallel")
#Set Up Model
weibull_cens_rand <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 2000,
                    sample = 5000, thin = 1, n.chains = 2)


weibull_cens_rand_mcmc <- as.mcmc.list(weibull_cens_rand)
# plot(weibull_cens_rand)

#####With Simulation of PPC####
#Bayes1_Weibull_random.r

# MODEL SPECIFICATION
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dweib(k, invlambda[i])
    #log(invlambda[i]) <- beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]]
    log(invlambda[i]) <- beta1 *x1[i] + b0[id[i]]
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dweib(k, invlambda[i])
    #log(invlambda[i]) <- beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]]
    log(invlambda[i]) <- beta1 *x1[i] + b0[id[i]]
  }
  #priors
  k ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  sigma ~ dunif(0,100)  
    
    for (i in 1:(N1+N2)){
    ppo[i] <- dweib(y[i],k, invlambda[i])
  }
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(mean(beta0 + beta2 *x2[i:(i+6)]) ,sigma)
    b0_rep[i] ~ dnorm(mean(beta0 + beta2 *x2[i:(i+6)]) ,sigma)
    #mu_gr[i] <- mean(beta0 + beta1 *x1[i] + beta2 *x2[i])
    #b0_rep[i]~ dnorm(0,sigma)
    
    
    # Ranked thetas
    rank_b01[i,1:Nsubj] <- sort(b0[1:Nsubj])
    rank_b0[i] <- rank_b01[i,i]
    rank_b0_rep1[i,1:Nsubj]  <- sort(b0_rep[1:Nsubj])
    rank_b0_rep[i] <- rank_b0_rep1[i,i]
  }
  # PPCs checking distribution of theta
  # min and max of residuals
  tmin1 <- sort(b0[])
  tmin <- tmin1[1]
  tmax1 <- sort(b0[])
  tmax <- tmax1[Nsubj]
  tmin.rep1 <- sort(b0_rep[])
  tmin.rep <-  tmin.rep1[1]
  tmax.rep1 <- sort(b0_rep[])
  tmax.rep <- tmax.rep1[Nsubj]
  tmin.test <- step(tmin.rep-tmin)
  tmax.test <- step(tmax.rep-tmax)
  # Kolmogorov-Smirnov test for residuals

  for (i in 1:Nsubj){
    F.gauss[i] <- phi( rank_b0[i])
    F.gauss.rep[i] <-  phi( rank_b0_rep[i])

    F.diff[i] <- max(F.gauss[i]-(i-1)/Nsubj,i/Nsubj-F.gauss[i])
    F.diff.rep[i] <- max(F.gauss.rep[i]-(i-1)/Nsubj,i/Nsubj-F.gauss.rep[i])
  }

  ks1 <- sort(F.diff[])
  ks <- ks1[Nsubj]
  ks.rep <- sort(F.diff.rep[])
  ks.rep2 <- ks.rep[Nsubj]

  ks.test <- step(ks.rep2-ks)

}"
parameters = c("ppo","k", "beta2", "beta0","sigma", "b0_rep", "beta1", "b0","tmin.test","tmax.test","ks.test", "sigma")


weibull_cens_rand_rep <- run.jags(model = model.function,
                         monitor = parameters, data = model.data,
                         inits = model.inits, burnin = 2000,
                         sample = 5000, thin = 1, n.chains = 2)

weibull_cens_rand_mcmc_rep <- as.mcmc.list(weibull_cens_rand_rep)


