setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library("runjags")
library(dplyr)
Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")
Grub <- Grub %>% arrange(upperlim)
length_Upper <- length(sort(Grub$upperlim))
#10
lenngth_NA_Upper <- nrow(Grub) - length_NA_Upper
NAs <- is.na(Grub$upperlim)
#just for numerical reasons, upper has to be bigger than lower
Grub$upperlim[NAs] <- 12.000001

#da alle Intervallcensored sind gerade muss ich immer eine 1 schicken
#https://stats.stackexchange.com/questions/13847/how-does-dinterval-for-interval-censored-data-work-in-jags
Grub$state <- c(rep(1,length_Upper),rep(2,lenngth_NA_Upper))
#also helpfull to understand later syntax: https://stats.stackexchange.com/questions/70858/right-censored-survival-fit-with-jags
Grub$value <- as.numeric(NA)
model.data <- list(y = Grub$value, z = Grub$state, N1 = length_Upper,N2 = lenngth_NA_Upper, x1 = Grub$grubsize,
                    x2 = Grub$group, lims = cbind(Grub$lowerlim,Grub$upperlim))


model.inits <- list(list(sigma=2, beta0=1, beta1 = 1,beta2 = 1 ),list(sigma=2, beta0=1, beta1 = 1,beta2 = 1 ))
parameters <-c("beta0", "beta1", "beta2", "sigma")
model.function <- "model{
  for (i in 1:N1){
    z[i] ~ dinterval(y[i], lims[i, ])
    y[i] ~ dlnorm(mu[i], sigma)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  for (i in (N1+1):(N1+N2)){
    z[i] ~ dinterval(y[i], lims[i,])
    y[i] ~ dlnorm(mu[i], sigma)
    mu[i] <- beta0 + beta1 *x1[i] + beta2 *x2[i]
  }
  #priors
  sigma ~ dunif(0.1,1000)
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
}"
runjags.options(method = "rjparallel")
ModelLogN <- run.jags(model = model.function,
                      monitor = parameters, data = model.data,
                      inits = model.inits, burnin = 5000, sample = 10000, thin = 1, n.chains = 2
                      , progress.bar = "text")


print(ModelLogN)



