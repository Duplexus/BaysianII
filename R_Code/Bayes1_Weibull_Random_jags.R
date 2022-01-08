#Bayes1_Weibull_random.r
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(124234)
source("helpfunctions.r")
library("runjags")
library("coda")
library("rjags")

Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")

model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group, id = Grub$id, Nsubj = length(unique(Grub$id)))
# model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
#                     id = Grub$id, Nsubj = length(unique(Grub$id)))

# MODEL SPECIFICATION 
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]]
  }
  #priors
  tau_b0 <- 1/sigma_b0
  sigma_b0 ~ dunif(0,100)
  k ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
     b0[i] ~ dnorm(0,tau_b0)
  }
}"

#Monitored Variables
parameters <- c("k", "sigma", "beta2", "beta0", "beta1","b0", "b0_rep","sigma_b0","y_sim")
#parameters <- c("k", "beta1", "b0", "b0_rep", "sigma")

# model.inits <- list(list(k=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20)),b0_rep = c(rep(1,times = 20))),
#                     list(k=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20)),b0_rep = c(rep(1,times = 20))))

model.inits <- list(list(k=4,sigma2_b0 =2, beta0=1, beta1 = 1,beta2 = 1,b0 = c(rep(1,times = 20)) ),
                    list(k=2,sigma2_b0 =20, beta0=10, beta1 = 10,beta2 = 10, b0 =  rnorm(20,0,30) ),
                    list(k=0.5,sigma2_b0 =15, beta0=20, beta1 = -10,beta2 = -15, b0 =  runif(20,0,10)  )
)

runjags.options(method = "rjparallel")
#Set Up Model
weibull_rand <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 2000,
                    sample = 5000, thin = 1, n.chains = 3)

weibull_rand_mcmc <- as.mcmc.list(weibull_rand)


#####With Simulation of PPC####
#Bayes1_Weibull_random.r

# MODEL SPECIFICATION
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dweib(k, invlambda[i])
    y_rep[i] ~ dweib(k, invlambda[i])
   invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]]
    #to generate the expected residuals
    y_sim[i] ~ dweib(k, invlambda[i])
    res[i] <- y[i] - y_sim[i]
  }
  #priors
  tau_b0 <- 1/sigma_b0
  sigma_b0 ~ dunif(0,100)
  #close to dgamma(1,0.000001) in bugs
  k ~ dunif(0,100)
  scale <- 1/k
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
    
    for (i in 1:N){
    ppo[i] <- dweib(y[i],k, invlambda[i])
    ppo_rep[i] <- dweib(y_rep[i],k, invlambda[i])
    D[i] <- -2*log(dweib(y[i],k, invlambda[i]))
    }
  Deviance <- sum(D[])
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
    b0_rep[i] ~ dnorm(0,tau_b0)
    
    
    # Ranked thetas also for PPC
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
parameters = c("ppo_rep","y_rep","res","ppo","k", "beta2", "beta0", "scale", "Deviance",
               "sigma", "b0_rep", "beta1", "b0","tmin.test","tmax.test","ks.test", "sigma","ss.test")


weibull_rand_rep <- run.jags(model = model.function,
                         monitor = parameters, data = model.data,
                         inits = model.inits, burnin = 2000,
                         sample = 5000, thin = 1, n.chains = 3)

weibull_rand_mcmc_rep <- as.mcmc.list(weibull_rand_rep)


