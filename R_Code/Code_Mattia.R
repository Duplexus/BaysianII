setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(coda)
library(rjags)
library(ggplot2)
library(R2WinBUGS)
data <- read.csv("..\\data\\Grubs_Nematodes.csv")

data
str(data)
summary(data)

data1 <- data
data1[4][is.na(data1[5])]
data1[5][is.na(data1[5])] <- 12

data1["MIDLIM"] <- apply(data1[c('LOWERLIM','UPPERLIM')], 1, mean)
data1

data1['Center_size'] = data1['GRUBSIZE'] - apply(data1['GRUBSIZE'], 2, mean)  
data1['Center_size']  = data1['Center_size'] / sd(as.vector(as.numeric(unlist(data1['Center_size']))))
#data1['Center_size'] <- data1['Center_size'] *100

#cbind(data1["Center_size"],data2,data2*sd(as.vector(as.numeric(unlist(data1['Center_size'])))))
#############################################
# 1 - Median death time per EPN
#############################################

group1 = data1$MIDLIM[data1['GROUP']==1]
group2 = data1$MIDLIM[data1['GROUP']==2]
median(group1) # =5
median(group2) # =3
data1$GROUP <- as.factor(data1$GROUP)

g1<-ggplot(data1, aes(x=MIDLIM, color=GROUP))
  g1 + 
    geom_histogram(binwidth = 1 ,fill='white', position = 'identity', alpha=0.4) +
    geom_vline(aes(xintercept=c(3),color='2')) +
    geom_vline(aes(xintercept=c(5),color='1'))
  
  
  
  
#############################################
# 2 - Effect of the covariate grubs'size
#############################################

#############################################
# Log-normal model
#############################################
  
set.seed(1234)

  model.data <- list(
    x = cbind(data1$GROUP,data1$Center_size), y = data1$MIDLIM, N = nrow(data1)
  )
  
  # MODEL SPECIFICATION 
  survival1 <- function(){
    # Specification data model
    for (i in 1:N)
    {
      y[i] ~ dlnorm(y.hat[i], tau)
      y.hat[i] <- beta0 + beta1*x[i,1]+beta2*x[i,2]
    }
    # Prior specification
    beta0 ~ dnorm(0.0, 0.000001)
    beta1 ~ dnorm(0.0, 0.000001)
    beta2 ~ dnorm(0.0, 0.000001)
    tau ~ dgamma(0.001, 0.001)
    sigma <- pow(tau, -2)
    # sigma ~ dunif(0,100)
  }
  write.model(survival1, "survival1.txt")
  
  # DEFINE INITIAL VALUES
  model.inits <- list(
    list(beta0 = 0, beta1 = 0, beta2 = 0, tau = 1 ),
    list(beta0 = -0.5, beta1 = -0.5, beta2 = -0.5, tau = 0.5),
    list(beta0 = 0.5, beta1 = 0.5, beta2 = 0.5, tau = 1.5)
  )
  
  # SET UP MODEL
  jags <- jags.model('survival1.txt',
                     data = model.data,
                     inits = model.inits,
                     n.chains = 3)
  
  update(jags, 10000)
  
  # Generate MCMC samples
  out <- coda.samples(jags,
                      c('beta0', 'beta1', 'beta2', 'tau'),
                      n.iter = 60000, thin=10)
  plot(out)
  
  autocorr.plot(out)
  
  # Posterior distributions
  densplot(out) 

  # Posterior summary statistics
  summary(out)
  
  dic1 <- dic.samples(jags, n.iter=60000, thin=10)
  dic1
  
#############################################
# Weibull AFT model
#############################################
  
  model.data <- list(
    x = cbind(data1$GROUP,data1$Center_size), y = data1$MIDLIM, N = nrow(data1)
  )
  
  # MODEL SPECIFICATION 
  survival2 <- function(){
    # Specification data model
    for (i in 1:N)
    {
      mu[i] <- beta0 + beta1*x[i,1]+beta2*x[i,2]
      inv_alpha[i] <- exp(-mu[i])
      lambda[i] <- pow(inv_alpha[i],gamma)
      y[i] ~ dweib(gamma,lambda[i])
    }
    # Prior specification
    beta0 ~ dnorm(0.0, 0.000001)
    beta1 ~ dnorm(0.0, 0.000001)
    beta2 ~ dnorm(0.0, 0.000001)
    gamma ~ dunif(0,10)
    sigma <- 1/gamma
  }
  write.model(survival2, "survival2.txt")
  
  # DEFINE INITIAL VALUES
  model.inits <- list(
    list(beta0 = 0, beta1 = 0, beta2 = 0, gamma = 2 ),
    list(beta0 = -0.5, beta1 = -0.5, beta2 = -0.5, gamma = 1.5),
    list(beta0 = 0.5, beta1 = 0.5, beta2 = 0.5, gamma = 2.5)
  )
  
  # SET UP MODEL
  jags <- jags.model('survival2.txt',
                     data = model.data,
                     inits = model.inits,
                     n.chains = 3)
  
  update(jags, 10000)
  
  # Generate MCMC samples
  out <- coda.samples(jags,
                      c('beta0', 'beta1', 'beta2', 'gamma'),
                      n.iter = 6000, thin=10)
  plot(out)
  
  autocorr.plot(out)
  
  # Posterior distributions
  densplot(out) 
  
  # Posterior summary statistics
  summary(out)
  
  dic2 <- dic.samples(jags, n.iter=60000, thin=10)
  
  # DIC log-normal : 636.7
  # DIC weibull : 641.3
  # I prefer the log-normal

  
#############################################
# 3 - Random effect
#############################################

####################################################
# Log-normal AFT model with random intercept plate 
####################################################
  
set.seed(420)
  
  model.data <- list(
    id = data1$UREPID, x = cbind(data1$GROUP,data1$Center_size), y = log(data1$MIDLIM),
    N = nrow(data1),  M = 20
  )
  
  # MODEL SPECIFICATION 
  survival3 <- function(){
    # Specification data model
    for (i in 1:N)
    {
      y[i] ~ dnorm(y.hat[i], tau)
      y.hat[i] <- b0[id[i]] + beta1*x[i,1]+beta2*x[i,2] + beta0
    }
    for(j in 1:M)
    {
      b0[j] ~ dnorm(0, tau_b0)
    }
    # Prior specification
    sigma ~ dunif(0,100)
    tau <- pow(sigma, -2)
    sigma_b0 ~ dunif(0,100)
    tau_b0 <- pow(sigma_b0, -2)
    
    beta0 ~ dnorm(0.0, 0.000001)
    beta1 ~ dnorm(0.0, 0.000001)
    beta2 ~ dnorm(0.0, 0.000001)
  }
  write.model(survival3, "survival3.txt")
  
  # DEFINE INITIAL VALUES
  model.inits <- list(
    list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = runif(1), sigma_b0 =  runif(1)),
    list(beta0 = -0.5, beta1 = -0.5, beta2 = -0.5, sigma = runif(1), sigma_b0 =  runif(1)),
    list(beta0 = 0.5, beta1 = 0.5, beta2 = 0.5, sigma = runif(1), sigma_b0 =  runif(1))
  )
  
  # SET UP MODEL
  jags <- jags.model('survival3.txt',
                     data = model.data,
                     inits = model.inits,
                     n.chains = 3)
  
  update(jags, 10000)
  
  # Generate MCMC samples
  out <- coda.samples(jags,
                      c('beta1', 'beta2', 'sigma', 'sigma_b0'),
                      n.iter = 6000, thin=10)
  plot(out)
  
  # Posterior distributions
  densplot(out) 
  
  # Posterior summary statistics
  summary(out)

  
#############################################
# 4 - Outlying observations
#############################################

  set.seed(420)
  
  model.data <- list(
    id = data1$UREPID, x = cbind(data1$GROUP,data1$Center_size), y = log(data1$MIDLIM),
    N = nrow(data1),  M = 20
  )
  
  # MODEL SPECIFICATION 
  survival3 <- function(){
    # Specification data model
    for (i in 1:N)
    {
      y[i] ~ dnorm(y.hat[i], tau)
      y.hat[i] <- beta0 + b0[id[i]] + beta1*x[i,1]+beta2*x[i,2]
      # CPO
      ppo[i] <- pow(2*3.141593,-0.5)*pow(sigma,-1)*exp(-0.5*pow((y[i]-y.hat[i])/sigma, 2))
      icpo[i] <- 1/ppo[i]
    }
    for(j in 1:M)
    {
      b0[j] ~ dnorm(0, tau_b0)
    }
    # Prior specification
    sigma ~ dunif(0,100)
    tau <- pow(sigma, -2)
    sigma_b0 ~ dunif(0,100)
    tau_b0 <- pow(sigma_b0, -2)
    
    beta0 ~ dnorm(0.0, 0.000001)
    beta1 ~ dnorm(0.0, 0.000001)
    beta2 ~ dnorm(0.0, 0.000001)
  }
  write.model(survival3, "survival3.txt")
  
  # DEFINE INITIAL VALUES
  model.inits <- list(
    list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = runif(1), sigma_b0 =  runif(1)),
    list(beta0 = -0.5, beta1 = -0.5, beta2 = -0.5, sigma = runif(1), sigma_b0 =  runif(1)),
    list(beta0 = 0.5, beta1 = 0.5, beta2 = 0.5, sigma = runif(1), sigma_b0 =  runif(1))
  )
  
  # SET UP MODEL
  jags <- jags.model('survival3.txt',
                     data = model.data,
                     inits = model.inits,
                     n.chains = 3)
  
  update(jags, 10000)
  
  # Generate MCMC samples
  out <- coda.samples(jags,
                      c('icpo'),
                      n.iter = 60000, thin=10)
  
  # Posterior summary statistics
  summary(out)

  icpo <- as.vector(summary(out)$statistics[,'Mean'])
  barplot(icpo, col="steelblue",ylim=c(0,80), ylab="1/CPO")
  
  which.max(icpo)
  which.max(icpo[-which.max(icpo)])
  icpo[64]
  icpo[57]
  geweke.plot(out,confidence = 0.95)
  # observations 57 and 64 are outliers
  
  
#############################################
# 5 - Distribution of the random effect
#############################################
  
  set.seed(420)
  data1['Center_size'] <- data1['Center_size'] *0.27
  model.data <- list(
    id = data1$UREPID, x = cbind(data1$GROUP,data1['Center_size']), y = log(data1$MIDLIM),
    N = nrow(data1),  M = 20
  )
  
  # MODEL SPECIFICATION 
  survival3 <- function(){
    # Specification data model
    for (i in 1:N)
    {
      y[i] ~ dnorm(y.hat[i], tau)
      y.hat[i] <- beta0 + b0[id[i]] + beta1*x[i,1]+beta2*x[i,2]
    }
    for(j in 1:M)
    {
      b0[j] ~ dnorm(0, tau_b0)
    }
    # Prior specification
    sigma ~ dunif(0,100)
    tau <- pow(sigma, -2)
    sigma_b0 ~ dunif(0,100)
    tau_b0 <- pow(sigma_b0, -2)
    
    beta0 ~ dnorm(0.0, 0.000001)
    beta1 ~ dnorm(0.0, 0.000001)
    beta2 ~ dnorm(0.0, 0.000001)
  }
  write.model(survival3, "survival3.txt")
  
  # DEFINE INITIAL VALUES
  model.inits <- list(
    list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = runif(1), sigma_b0 =  runif(1)),
    list(beta0 = -0.5, beta1 = -0.5, beta2 = -0.5, sigma = runif(1), sigma_b0 =  runif(1)),
    list(beta0 = 0.5, beta1 = 0.5, beta2 = 0.5, sigma = runif(1), sigma_b0 =  runif(1))
  )
  
  # SET UP MODEL
  jags <- jags.model('survival3.txt',
                     data = model.data,
                     inits = model.inits,
                     n.chains = 3)
  
  update(jags, 10000)
  
  # Generate MCMC samples
  out <- coda.samples(jags,
                      c('b0',"beta0","beta1","beta2"),
                      n.iter = 60000, thin=10)
  
  
  # Posterior summary statistics
  summary(out)
  
  # histogram and qqplot of posterior means of the random intercepts
  
  postmean_b0 <- as.vector(summary(out)$statistics[,'Mean'])
  hist(postmean_b0)
  qqnorm(postmean_b0)
  qqline(postmean_b0)
  
  
  # histogram and qqplot of the random intercepts
  
  b0 <- as.vector(c(as.vector(out[[1]]), as.vector(out[[2]]), as.vector(out[[3]])))
  hist(b0, freq=F)
  lines(seq(-3,2,0.001),dnorm(seq(-3,2,0.001),0,sd(b0)), col=2)
  qqnorm(b0)
  qqline(b0)
  
  
