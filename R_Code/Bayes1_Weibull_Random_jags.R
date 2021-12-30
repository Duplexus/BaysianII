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

# MODEL SPECIFICATION 
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dweib(k, invlambda[i])
    log(invlambda[i]) <- beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]]
  }
  #priors
  k ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,0.1)
  }
}"

#Monitored Variables
parameters <- c("k", "beta2", "beta0", "beta1", "b0")

model.inits <- list(list(k=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20))),
                    list(k=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20))))


runjags.options(method = "rjparallel")
#Set Up Model
weibull_rand <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 2000,
                    sample = 5000, thin = 1, n.chains = 2)

# 
# Weibull_bayes <- read.bugs(model.out)
# 
# summary_weibull <- summary(Weibull_bayes)
# b_params <- summary_weibull$statistics[grepl("b0\\[", dimnames(summary_weibull$statistics)[[1]]),]
# mean(b_params[1:10,1])
# attributes(summary_weibull$statistics)
# 
# 
#####With Simulation of PPC####
#Bayes1_Weibull_random.r

Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")

model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group, id = Grub$id, Nsubj = length(unique(Grub$id)))

# MODEL SPECIFICATION
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dweib(k, invlambda[i])
    log(invlambda[i]) <- beta0 + beta1 *x1[i] + beta2 *x2[i] + b0[id[i]]
  }
  #priors
  k ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
    for (i in 1:N){
    ppo[i] <- dweib(y[i],k, invlambda[i])
  }
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,0.1)
    b0_rep[i]~ dnorm(0,0.1)
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
parameters = c("ppo","k", "beta2", "beta0", "beta1", "b0","tmin.test","tmax.test","ks.test")


weibull_rand_rep <- run.jags(model = model.function,
                         monitor = parameters, data = model.data,
                         inits = model.inits, burnin = 2000,
                         sample = 5000, thin = 1, n.chains = 2)

weibull_rand_mcmc_rep <- as.mcmc.list(weibull_rand_rep)

# 
# Weibull_bayes <- read.bugs(model.out)
# 
# summary_weibull <- summary(Weibull_bayes)
# b_params <- summary_weibull$statistics[grepl("b0\\[", dimnames(summary_weibull$statistics)[[1]]),]
# mean(b_params[1:10,1])
# attributes(summary_weibull$statistics)
# #how unnormal the random effects
# mean(get_values(Weibull_bayes,"tmax.test"))
# mean(get_values(Weibull_bayes,"tmin.test"))
# mean(get_values(Weibull_bayes,"ks.test"))
# 



