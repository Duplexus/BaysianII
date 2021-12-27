#Bayes1_Weibull_random.r
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
set.seed(124234)
library(coda)
library(R2OpenBUGS)

Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")

model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group, id = Grub$id, Nsubj = length(unique(Grub$id)))

# MODEL SPECIFICATION 
model.function <- function(){
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
}
write.model(model.function, "Scripts\\Bayes1_Weibull_Random.txt")
model.inits <- function(){list(k=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20)) )}
parameters = c("k", "beta2", "beta0", "beta1", "b0")


model.out <- bugs(model.data, model.inits, 
                  model.file = "Bayes1_Weibull_Random.txt",
                  parameters=parameters,
                  n.chains = 2, n.iter = 5000,n.thin = 100,  n.burnin = 2000, debug = T,
                  codaPkg=T,
                  working.directory = ".\\Scripts")

Weibull_bayes <- read.bugs(model.out)

summary_weibull <- summary(Weibull_bayes)
b_params <- summary_weibull$statistics[grepl("b0\\[", dimnames(summary_weibull$statistics)[[1]]),]
mean(b_params[1:10,1])
attributes(summary_weibull$statistics)

#####With Simulation of PPC####
#Bayes1_Weibull_random.r

Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")

model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group, id = Grub$id, Nsubj = length(unique(Grub$id)))

# MODEL SPECIFICATION 
model.function <- function(){
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
    b0_rep[i]~ dnorm(0,0.1)
    # Ranked thetas
    rank_b0[i] <- ranked( b0[1:Nsubj],i )
    rank_b0_rep[i] <- ranked( b0_rep[1:Nsubj],i ) 
  }
  # PPCs checking distribution of theta
  # min and max of residuals
  tmin <- ranked(b0[],1)
  tmax <- ranked(b0[],Nsubj)
  tmin.rep <- ranked(b0_rep[],1)
  tmax.rep <- ranked(b0_rep[],Nsubj)
  tmin.test <- step(tmin.rep-tmin)
  tmax.test <- step(tmax.rep-tmax)
  # Kolmogorov-Smirnov test for residuals
  
  for (i in 1:Nsubj){
    F.gauss[i] <- phi( rank_b0[i])
    F.gauss.rep[i] <-  phi( rank_b0_rep[i])
    
    F.diff[i] <- max(F.gauss[i]-(i-1)/Nsubj,i/Nsubj-F.gauss[i])
    F.diff.rep[i] <- max(F.gauss.rep[i]-(i-1)/Nsubj,i/Nsubj-F.gauss.rep[i])
  }  
  
  ks <- ranked(F.diff[],Nsubj)
  ks.rep <- ranked(F.diff.rep[],Nsubj)
  
  ks.test <- step(ks.rep-ks)  
  
}
write.model(model.function, "Scripts\\Bayes1_Weibull_Random.txt")
model.inits <- function(){list(k=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20)) )}
parameters = c("k", "beta2", "beta0", "beta1", "b0","tmin.test","tmax.test","ks.test")


model.out <- bugs(model.data, model.inits, 
                  model.file = "Bayes1_Weibull_Random.txt",
                  parameters=parameters,
                  n.chains = 2, n.iter = 5000,n.thin = 1,  n.burnin = 2000, debug = T,
                  codaPkg=T,
                  working.directory = ".\\Scripts")

Weibull_bayes <- read.bugs(model.out)

summary_weibull <- summary(Weibull_bayes)
b_params <- summary_weibull$statistics[grepl("b0\\[", dimnames(summary_weibull$statistics)[[1]]),]
mean(b_params[1:10,1])
attributes(summary_weibull$statistics)
#how unnormal the random effects
mean(get_values(Weibull_bayes,"tmax.test"))
mean(get_values(Weibull_bayes,"tmin.test"))
mean(get_values(Weibull_bayes,"ks.test"))




