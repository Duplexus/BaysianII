setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("helpfunctions.r")
library("runjags")
library("coda")
library("rjags")
Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")
#attributes(Grub$grubsize) <- NULL


#DEFINE INTITIAL VALUES
model.inits <- list(list(tau=2, beta0=1, beta1 = 1,beta2 = 1 ),
                    list(tau=20, beta0=10, beta1 = 10,beta2 = 10 ),
                    list(tau=15, beta0=20, beta1 = -10,beta2 = -15 )
                    )

#Monitored Variables
parameters <-c("beta0", "beta1", "beta2", "sigma","tau")



#what happens if big grub size is forbidden



#Initial Values
model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group)

# Specification data model
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  #priors
  tau ~ dgamma(0.001, 0.001)
  sigma <- sqrt(1/tau)
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
}"

runjags.options(method = "rjparallel")
#Set Up Model
#Generate MCMC SAMpls
lognorm <- run.jags(model = model.function,
                     monitor = parameters, data = model.data,
                     inits = model.inits, burnin = 2000,
                     sample = 5000, thin = 1, n.chains = 3)
lognorm_mcmc <- as.mcmc.list(lognorm)
    # ####If not used in Master decomment once####
    # #Model Diagnostics and so on 
    # #Starts with ModelLogN and plots it and calculates BIC and the 
    # #Statistics of Rubin and so on
    # #some Results
    # plot(lognorm)
    # print(lognorm)
    # #The DIC Value for model comparison
    # dic_val <- extract.runjags(lognorm, "dic")
    # dic_val
    # extract.runjags(ModelLogN, "stochastic")
    # 
    # #coda integration so also coda stuff is available
    # #Model Diagnostik plots
    # mcmc_lognorm <- as.mcmc.list(lognorm)
    # 
    # gelman.diag(mcmc_lognorm, confidence = 0.95)
    # gelman.plot(mcmc_lognorm, confidence = 0.95)
    # geweke.diag(mcmc_lognorm)
    # geweke.plot(mcmc_lognorm)
    # plot(mcmc_lognorm)


####Model Disagnostics and PPC####
model.function <- "model{
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
    y_rep[i] ~ dlnorm(mu[i], tau)
  }
  #priors
  tau ~ dgamma(0.001, 0.001)
  sigma <- sqrt(1/tau)
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
  #ppc look on residuals
  for (i in 1:N){
    ppo[i] <- dlnorm(y[i],mu[i],tau)
    ppo_rep[i] <- dlnorm(y_rep[i],mu[i],tau)
    k[i] <- log(y[i])
    res[i] <- k[i] - mu[i]
    #for DIC
    D[i] <- -2*log(dlnorm(y[i],mu[i], tau))
  }
  Deviance <- sum(D[])
}"
parameters <-c("ppo_rep","beta0", "beta1", "beta2", "sigma","ppo","res","mu",
               "Deviance","y_rep","tau")
lognormal_rep <- run.jags(model = model.function,
                      monitor = parameters, data = model.data,
                      inits = model.inits, burnin = 2000,
                      sample = 5000, thin = 10, n.chains = 3)
lognorm_mcmc_rep <- as.mcmc.list(lognormal_rep)



# #1.) extract the ppo values 
# subset_pred <- grepl("ppo\\[", dimnames(lognormal_mcmc_rep[[1]])[[2]])
# mcmc_subset <- get_values(lognormal_mcmc_rep,subset_pred)
# #ppo which are now far of?
# plot(1/apply(as.matrix(mcmc_subset),2,mean))
# 
# #2.) extract the res values 
# subset_pred <- grepl("res\\[", dimnames(lognormal_mcmc_rep[[1]])[[2]])
# mcmc_subset <- get_values(lognormal_mcmc_rep,subset_pred)
# #how look the average residuals in the log world?
# hist(apply(mcmc_subset,2,mean), breaks = 20)
# 
# #3.) Get the DIC running
# subset_pred <- grepl("Deviance", dimnames(lognormal_mcmc_rep[[1]])[[2]])
# mcmc_subset <- get_values(lognormal_mcmc_rep,subset_pred)
# mcmc_subset_dic<- mcmc_subset
# md <- mean(mcmc_subset_dic)
# 
# a <- summary(lognormal_rep)
# a <- a[c("beta0","beta1","beta2","sigma"),"Mean"]
# pd <- md - (-2 *sum(log(dlnorm(Grub$value,a["beta0"]+a["beta1"]*Grub$grubsize 
#                    + a["beta2"]*Grub$group,sqrt(1/a["sigma"])))))
# c(pd,md,pd+md)
# dic_val <- extract.runjags(lognormal_rep, "dic")
# dic_val
# summary(dic_val)



# library(ggplot2)
# summary(lm(I(log(Grub$value)) ~ Grub$grubsize * Grub$group))
# Grub$group <- as.factor()
# ggplot(Grub, aes(y = value, x = grubsize, colour = as.factor(group))) + geom_point()
# 


