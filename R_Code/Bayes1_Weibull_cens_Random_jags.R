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

parameters = c("sigma", "beta2", "beta0", "beta1", "b0","sigma2_b0")

model.data <- list(y = Grub$value2, z = Grub$state, N1 = length_Upper,N2 = lenngth_NA_Upper,
                   x1 = Grub$grubsize, x2 = Grub$group, id = Grub$id,
                   lims = cbind(Grub$lowerlim,Grub$upperlim),
                   Nsubj = length(unique(Grub$id)))
# MODEL SPECIFICATION 
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dweib(k, invlambda[i])
    log(invlambda[i]) <- -(beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]])
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]]
  }
  #priors
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  # similar to dgamma(1,0.0001) bzw. dgamma(1,10000) if interpreted as we would get it from the exponential 
  k ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
  }
}"

#Monitored Variables
parameters <- c("k", "scale", "beta2", "beta0", "beta1","b0","sigma2_b0")
#parameters <- c("k", "beta1", "b0", "b0_rep", "sigma")

# model.inits <- list(list(k=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20)),b0_rep = c(rep(1,times = 20))),
#                     list(k=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20)),b0_rep = c(rep(1,times = 20))))
model.inits <- list(list(k=4,sigma2_b0 = 2, beta0=1, beta1 = 1,beta2 = 1,b0 = c(rep(1,times = 20)) ),
                    list(k=2,sigma2_b0 = 20, beta0=10, beta1 = 10,beta2 = 10, b0 =  rnorm(20,0,4) ),
                    list(k=0.5,sigma2_b0 = 15, beta0=20, beta1 = -10,beta2 = -15, b0 =  runif(20,-5,5)  )
)

runjags.options(method = "rjparallel")
#Set Up Model
weibull_cens_rand <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 10000,
                    sample = 5000, thin = 10, n.chains = 3)


weibull_cens_rand_mcmc <- as.mcmc.list(weibull_cens_rand)
# plot(weibull_cens_rand)

#####With Simulation of PPC####
#Bayes1_Weibull_random.r

# MODEL SPECIFICATION
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dweib(k, invlambda[i])
    y_rep[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]]
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dweib(k, invlambda[i])
    y_rep[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]]
    
  }
  #priors
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  k ~ dunif(0,100)
  scale <- 1/k
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  sigma ~ dunif(0,100)  
    
    for (i in 1:(N1+N2)){
    ppo[i] <- dweib(y[i],k, invlambda[i])
    ppo_rep[i] <- dweib(y_rep[i],k, invlambda[i])
    #for DIC
    D[i] <- -2*log(dweib(y[i],k, invlambda[i]))
    }
  Deviance <- sum(D[])
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
    b0_rep[i] ~ dnorm(0,tau_b0)
    
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
  
  # Sinharay and Stern test
  #hard coding since we have 20 random effects this is fine
  nmed <- round(20/2)
  tmed1 <- sort(b0[])
  tmed <- tmed1[nmed]
  tmed.rep1 <- sort(b0_rep[])
  tmed.rep <- tmed1[nmed]

  ss <- abs(tmax-tmed)-abs(tmin-tmed)
  ss.rep <-abs(tmax.rep-tmed.rep)-abs(tmin.rep-tmed.rep)
  ss.test <- step(ss.rep-ss) 

}"
parameters = c("y_rep","ppo_rep","ppo","k", "beta2", "beta0","sigma", "b0_rep", "beta1", "b0",
               "tmin.test","tmax.test","ks.test", "sigma","Deviance","scale","ss.test","y","sigma2_b0")


weibull_cens_rand_rep <- run.jags(model = model.function,
                         monitor = parameters, data = model.data,
                         inits = model.inits, burnin = 2000,
                         sample = 5000, thin = 1, n.chains = 3)

weibull_cens_rand_mcmc_rep <- as.mcmc.list(weibull_cens_rand_rep)


