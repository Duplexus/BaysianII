setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library("runjags")
Grub <- read.csv("..\\data\\Grubs_Easy.csv")
NAs <- is.na(Grub$upperlim)
Grub$upperlim[NAs]<- 14
Grub$value <- 1
Grub$value  <- as.numeric(NA)
model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group, lims = cbind(Grub$lowerlim,Grub$upperlim))



model.inits <- list(sigma=2, beta0=1, beta1 = 1,beta2 = 1 )
parameters <-c("beta0", "beta1", "beta2", "sigma")
model.function <- "model{
  for (i in 1:N){
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
                      inits = model.inits, burnin = 2000, sample = 5000, thin = 1, n.chains = 2)


print(ModelLogN)



