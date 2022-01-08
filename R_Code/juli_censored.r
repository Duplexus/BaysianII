setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(coda)
library(rjags)
library(ggplot2)
library(R2WinBUGS)
library(runjags)

data <- read.csv("..\\data\\Grubs_Nematodes.csv")
data
time <- data$LOWERLIM
cens <- matrix(c(time, rep(NA, length(time))), nrow = length(time), ncol = 2)
cens[,1] <- data$LOWERLIM
cens[,2] <- data$UPPERLIM
cens

right_censored <- which(cens[,2] %in% NA)
interval_censored <- 1:nrow(data) 
interval_censored <- interval_censored[-right_censored]

is.censored <- rep(1, length(data$UPPERLIM))
is.censored  

time <- rep(NA, length(data$UPPERLIM))
time



model.data <- list(
  x = cbind(data$GROUP, scale(data$GRUBSIZE)), y = time, right_censored = right_censored,
  interval_censored = interval_censored, is.censored = is.censored, cens = cens
)




survival1 <- function(){
  # Specification data model
  for (i in interval_censored)
  {
    is.censored[i] ~ dinterval(y[i], cens[i,])
    y[i] ~ dlnorm(y.hat[i], tau)
    y.hat[i] <- beta0 + beta1*x[i,1]+beta2*x[i,2]
  }
  for (i in right_censored)
  {
    is.censored[i] ~ dinterval(y[i], cens[i,1])
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

update(jags, 1000)


# Generate MCMC samples
out_censored <- coda.samples(jags,
                             c('beta0', 'beta1', 'beta2', 'tau'),
                             n.iter = 4000, thin=10)
plot(out_censored)

# Posterior distributions
densplot(out) 

# Posterior summary statistics
summary(out_censored)

dic1 <- dic.samples(jags, n.iter=5000, thin=10)
dic1
