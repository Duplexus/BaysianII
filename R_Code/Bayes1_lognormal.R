#lognormal Modell 1 no random effects
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(coda)
library(R2OpenBUGS)
Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")
# Grub$grubsize <- scale(Grub$grubsize)
# attributes(Grub$grubsize) <- NULL
# # Grub$group <- scale(Grub$group)

model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group)

# MODEL SPECIFICATION 
model.function <- function(){
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], sigma)
    predict[i] ~ dlnorm(mu[i], sigma)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
    
  
  }
  #priors
  sigma ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
}
write.model(model.function, "Scripts\\Bayes1_lognormal.txt")
model.inits <- function(){list(sigma=2, beta0=1, beta1 = 1,beta2 = 1, predict = c(rep(0,times = 140)))}
# parameters = c("sigma", "beta2", "beta0", "beta1")

parameters = c("sigma", "beta2", "beta0", "beta1","predict")

model.out <- bugs(model.data, model.inits, 
                  model.file = "Bayes1_lognormal.txt",
                  parameters=parameters,
                  n.chains = 2, n.iter = 5000,  n.burnin = 0, debug = T,
                  codaPkg=T,
                  working.directory = ".\\Scripts")
#####Results####

lognormal1 <- read.bugs(model.out)
summary(lognormal1)

#parametric
library(survival)
lognormal_results1 <- survreg(Surv(value) ~ grubsize +  group, Grub, dist = "lognormal")
# scale is sigma seems to be 1/the one from the baysian but beside this similar
#results between bayes and model are really similar
summary(lognormal_results1)

## Model is runned with predictions again
####PPC####
#wie sieht das geschätzte Modell aus:
subset_pred <- grepl("predict\\[", dimnames(lognormal1[[1]])[[2]])
#extracts all values and packages them into one list
#... one can include which rows are wanted
#object is object after coda
get_values <- function(object,column,...){
  for ( i in 1:length(object)){
    object[[i]] <- object[[i]][,column]
  }
  return(as.matrix(object)[...,])
}

#applys a function on every row and compares it to the original value
#the function always has to reduce the dim from n -> 1
vgl_fct <- function(x,fct,orig){
  result <- apply(x,1,fct)
  return(sum(result > orig)/length(result))
}


data_conv <- (get_values(lognormal1,subset_pred))
#mean similar?
vgl_fct(data_conv,fct = mean, mean(Grub$value))
#sd similar
vgl_fct(data_conv,fct = sd, sd(Grub$value))

#wie viel Prozent der Residuuen sind über zwei Standardeinheiten weg?
#takes a vector
outer_5 <- function(x,sds=1.96){
  l <- length(x)
  a1 <- x >= sds*sd(x) + mean(x)
  a2 <- x <= - sds*sd(x) + mean(x)
  return(list(left = sum(a2)/l, right = sum(a1)/l, comb = sum(a2)/l+ sum(a1)/l))

}
#just log because we live in the log world
data <- log(data_conv[200,]) - log(Grub$value)
outer_5(data)$comb
hist(data)

















#lognormal Modell 1 no random effects
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(coda)
library(R2OpenBUGS)
Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")
# Grub$grubsize <- scale(Grub$grubsize)
# attributes(Grub$grubsize) <- NULL
# # Grub$group <- scale(Grub$group)

model.data <- list( y = (Grub$value), N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group)

# MODEL SPECIFICATION 
model.function <- function(){
  for (i in 1:N){
    y[i] ~ dnorm(mu[i], sigma)
    predict[i] ~ dnorm(mu[i], sigma)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
    D[i] <- - log(sigma) + log(2*3.14159265358979) + pow(log(y[i])-mu[i],2)*sigma
    
  }
  #priors
  sigma ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  Deviance <- sum(D[])
}
write.model(model.function, "Scripts\\Bayes1_lognormal.txt")
model.inits <- function(){list(sigma=2, beta0=1, beta1 = 1,beta2 = 1, predict = c(rep(0,times = 140)))}
# parameters = c("sigma", "beta2", "beta0", "beta1")

parameters = c("sigma", "beta2", "beta0", "beta1","predict","Deviance")

model.out <- bugs(model.data, model.inits, 
                  model.file = "Bayes1_lognormal.txt",
                  parameters=parameters,
                  n.chains = 2, n.iter = 5000,  n.burnin = 0, debug = T,
                  codaPkg=T,
                  working.directory = ".\\Scripts")
#####Results####


#3.) Get the DIC running
subset_pred <- grepl("Deviance", dimnames(lognormal1[[1]])[[2]])
mcmc_subset <- get_values(lognormal1,subset_pred)
mean(mcmc_subset)
