#Weibull Modell 1 no random effects
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(coda)
library(R2OpenBUGS)
Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")

model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group)


model.inits <- list(list(k=3, beta0=1, beta1 = 1,beta2 = 1 ),
                    list(k=1, beta0=10, beta1 = 10,beta2 = 10 ),
                    list(k=0.5, beta0=20, beta1 = -10,beta2 = -15 )
)


# MODEL SPECIFICATION 
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
    # Distribution of future observed counts for plate i
    #predict[i]  ~ dweib(k, invlambda[i]) 
  }
  #priors
  scale <- 1/k
  #similar to a uniform prior, but goes till ifnitite
  k ~ dgamma(1,0.01)
  beta0 ~ dnorm(0,10e-1)
  beta1 ~ dnorm(0,10e-1)
  beta2 ~ dnorm(0,10e-1)
}"


#Monitored Variables
parameters <-c("beta0", "beta1", "beta2","scale")

runjags.options(method = "rjparallel")
#Set Up Model
weibull_sens <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 2000,
                    sample = 5000, thin = 1, n.chains = 3)
# MODEL SPECIFICATION 
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
    # Distribution of future observed counts for plate i
    #predict[i]  ~ dweib(k, invlambda[i]) 
  }
  #priors
  scale <- 1/k
  #really high weight on really small values and then fast decreasing
  k ~ dgamma(1,1)
  beta0 ~ dnorm(0,10e-5)
  beta1 ~ dnorm(0,10e-5)
  beta2 ~ dnorm(0,10e-5)
}"


#Monitored Variables
parameters <-c("beta0", "beta1", "beta2","scale")

runjags.options(method = "rjparallel")
#Set Up Model
weibull_sens2 <- run.jags(model = model.function,
                         monitor = parameters, data = model.data,
                         inits = model.inits, burnin = 2000,
                         sample = 5000, thin = 1, n.chains = 3)

