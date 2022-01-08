#Weibull Modell 1 no random effects
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(coda)
library(R2OpenBUGS)
Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")

model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group)

# MODEL SPECIFICATION 
#now the model specification is in line with the estimation in survival
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
#Set Up Model
weibull <- run.jags(model = model.function,
                    monitor = parameters, data = model.data,
                    inits = model.inits, burnin = 5000,
                    sample = 10000, thin = 1, n.chains = 3)

# weibull_2 <- as.mcmc.list(weibull)
# #k = 1.6
# Weibull_summary <- summary(weibull_2)
# Weibull_summary$statistics
weibull_mcmc <- as.mcmc.list(weibull)



model.function <- "model{
  for (i in 1:N){
    y[i] ~ dweib(k, invlambda[i])
    y_rep[i] ~ dweib(k, invlambda[i])
    invlambda[i] <- pow(t[i], k)
    t[i] <- exp(-h[i])
    h[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
    # Distribution of future observed counts for plate i
    #predict[i]  ~ dweib(k, invlambda[i]) 
  }
  #priors
  scale <- 1/k
  k ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  for (i in 1:N){
    ppo[i] <- dweib(y[i],k, invlambda[i])
    #for DIC
    D[i] <- -2*log(dweib(y[i],k, invlambda[i]))
  }
  Deviance <- sum(D[])
}"

parameters <-c("beta0", "beta1", "beta2", "invlambda","scale","ppo","Deviance"
               ,"y_rep")
weibull_rep <- run.jags(model = model.function,
                          monitor = parameters, data = model.data,
                          inits = model.inits, burnin = 5000,
                          sample = 10000, thin = 1, n.chains = 3)

weibull_mcmc_rep <- as.mcmc.list(weibull_rep)






# 
# 
# #parametric 
# library(survival)
# estimates <- survreg(Surv(value) ~ grubsize +  group, Grub, dist = "weibull")
# estimates
# 1/0.62 # same in the weibull case if you flip the sign
# Weibull_summary$statistics
