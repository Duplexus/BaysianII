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


model.inits <- list(list(tau=2,sigma2_b0 =2 , beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20))),
                    list(tau=20,sigma2_b0 =20, beta0=10, beta1 = 10,beta2 = 10 , b0 =  rnorm(20,0,30) ),
                    list(tau=15,sigma2_b0 =15, beta0=20, beta1 = -10,beta2 = -15 , b0 =  runif(20,0,10))
)
parameters = c("sigma", "beta0", "beta1", "beta2","sigma2_b0","b0")

model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    id = Grub$id,x2 = Grub$group,  Nsubj = length(unique(Grub$id)))



model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    
    
  }
  #priors
  tau ~ dgamma(0.001,0.001)
  sigma <- sqrt(1/tau)
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

lognorm_rand <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 2000,
                    sample = 5000, thin = 1, n.chains = 3)
lognorm_rand_mcmc <- as.mcmc.list(lognorm_rand)

summary(lognorm_rand)
summary(lognorm_rand_mcmc)
#parameters = c("sigma", "beta2", "beta0", "beta1", "b0","predict","b0.rep","ppo")
parameters = c("sigma", "beta0","beta1", "beta2", "b0","predict","b0.rep","ppo","sigma2_b0","Deviance",
               "tmax.test","tmin.test","ks.test","ss.test","y_rep", "ppo_rep","ppo2","tau")
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    y_rep[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]+ b0[id[i]]
    #for DIC
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
    ppo[i] <- dlnorm(y[i],mu[i],tau)
   #Explicitly ppo
    ppo2[i] <- pow(tau/(2*3.141593),0.5)*pow(y[i],-1)*exp(-tau*(pow((log(y[i])-mu[i]), 2))/2)

  }
  Deviance <- sum(D[])
  #priors
  tau ~ dgamma(0.001, 0.001)
  #now it is from inverse gamma
  sigma <- sqrt(pow(tau,-1))
  
  tau_b0 <- 1/sigma2_b0
  sigma2_b0 ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  for (i in 1:N){
    ppo_rep[i] <- dlnorm(y_rep[i],mu[i],tau)
  }
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,tau_b0)
    # Distribution of future b0_i
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


lognorm_rand_rep <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 2000,
                    sample = 5000, thin = 1, n.chains = 3)
lognorm_rand_mcmc_rep <- as.mcmc.list(lognorm_rand_rep)

