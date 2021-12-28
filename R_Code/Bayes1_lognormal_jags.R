setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("helpfunctions.r")
library("runjags")
library("coda")
library("rjags")
Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")
#attributes(Grub$grubsize) <- NULL

#DEFINE INTITIAL VALUES
model.inits <- list(list(sigma=2, beta0=1, beta1 = 1,beta2 = 1 ),
                    list(sigma=2, beta0=1, beta1 = 1,beta2 = 1 )
                    )

#Monitored Variables
parameters <-c("beta0", "beta1", "beta2", "sigma")

#Initial Values
model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group)

# Specification data model
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], sigma)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  #priors
  sigma ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
}"

runjags.options(method = "rjparallel")
#Set Up Model
#Generate MCMC SAMpls
ModelLogN <- run.jags(model = model.function,
                     monitor = parameters, data = model.data,
                     inits = model.inits, burnin = 2000,
                     sample = 5000, thin = 1, n.chains = 2)

#results <- extend.jags(ModelLogN, sample=5000)



#Model Diagnostics and so on 
#Starts with ModelLogN and plots it and calculates BIC and the 
#Statistics of Rubin and so on
#some Results
plot(ModelLogN)
print(ModelLogN)
#The DIC Value for model comparison
dic_val <- extract.runjags(ModelLogN, "dic")
dic_val
extract.runjags(ModelLogN, "stochastic")

#coda integration so also coda stuff is available
#Model Diagnostik plots
mcmc <- as.mcmc.list(ModelLogN)

gelman.diag(mcmc, confidence = 0.95)
gelman.plot(mcmc, confidence = 0.95)
geweke.diag(mcmc)
geweke.plot(mcmc)
plot(mcmc)


#PPC for normality of Errors
####Model Disagnostics and PPC
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dnorm(mu[i], sigma)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
    ppo[i] <- dnorm(k[i],mu[i],sigma)
    k[i] <- log(y[i])
  }
  #priors
  sigma ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
  #ppc look on residuals
  for (i in 1:N){
    res[i] <- k[i] - mu[i]
    D[i] <- dnorm(y[i],mu[i], sigma)
  }
  Deviance <- sum(D[])
}"
parameters <-c("beta0", "beta1", "beta2", "sigma","ppo","res","mu","Deviance")
ModelLogN_rep <- run.jags(model = model.function,
                      monitor = parameters, data = model.data,
                      inits = model.inits, burnin = 2000,
                      sample = 5000, thin = 1, n.chains = 2)
mcmc_rep <- as.mcmc.list(ModelLogN_rep)


#1.) extract the ppo values 
subset_pred <- grepl("ppo\\[", dimnames(mcmc_rep[[1]])[[2]])
mcmc_subset <- get_values(mcmc_rep,subset_pred)
#ppo which are now far of?
plot(1/apply(as.matrix(mcmc_subset),2,mean))

#2.) extract the res values 
subset_pred <- grepl("res\\[", dimnames(mcmc_rep[[1]])[[2]])
mcmc_subset <- get_values(mcmc_rep,subset_pred)
#how look the average residuals in the log world?
hist(apply(mcmc_subset,2,mean), breaks = 20)

#3.) Get the DIC running
subset_pred <- grepl("D\\[", dimnames(mcmc_rep[[1]])[[2]])
mcmc_subset <- get_values(mcmc_rep,subset_pred)
erstes <- cbind(apply(mcmc_subset,2,mean))


subset_pred <- grepl("Deviance", dimnames(mcmc_rep[[1]])[[2]])
mcmc_subset <- get_values(mcmc_rep,subset_pred)
-2*mean(mcmc_subset)

dic_val <- extract.runjags(ModelLogN_rep, "dic")
dic_val











subset_pred <- grepl("mu\\[", dimnames(mcmc_rep[[1]])[[2]])
mcmc_subset_mean <- get_values(mcmc_rep,subset_pred)
subset_pred <- grepl("sigm", dimnames(mcmc_rep[[1]])[[2]])
mcmc_subset_sd <- as.vector(get_values(mcmc_rep,subset_pred))


subset_pred <- grepl("ppo\\[", dimnames(mcmc_rep[[1]])[[2]])
mcmc_subset <- get_values(mcmc_rep,subset_pred)

k <- -2 *apply(log(mcmc_subset),2,mean)
0.5 * var(k) +mean(k)

sum(log((as.matrix(mcmc_subset))))/nrow(mcmc_subset)

cbind(dic_val[["deviance"]],-2 *apply(log(mcmc_subset),2,mean),erstes)

dic_val <- extract.runjags(ModelLogN_rep, "dic")
dic_val
