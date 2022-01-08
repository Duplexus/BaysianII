#Bayes1_lognormal_random.r
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(coda)
library(R2OpenBUGS)
library(tidyverse)


"\\ /"
Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")
Grub <- Grub %>% arrange(upperlim)
length_Upper <- length(sort(Grub$upperlim))
#10
lenngth_NA_Upper <- nrow(Grub) - length_Upper
NAs <- is.na(Grub$upperlim)
#just for numerical reasons, second has to be bigger than first
#and just to lazy to separate therefore makes no difference
Grub$upperlim[NAs] <- 12.000001
Grub$state <- c(rep(1,length_Upper),rep(2,lenngth_NA_Upper))
Grub$value2 <- as.numeric(NA)


Grub$value

# model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
#                     x2 = Grub$group, id = Grub$id, Nsubj = length(unique(Grub$id)))

model.data <- list(N1 = length_Upper,N2 = lenngth_NA_Upper,
                   x1 = Grub$grubsize, x2 = Grub$group, id = Grub$id,
                   lims = cbind(Grub$lowerlim,Grub$upperlim),
                   Nsubj = length(unique(Grub$id) ))


# MODEL SPECIFICATION 





model.function <- function(){
  for (i in 1:N1){
    mu[i] <- beta0 + beta1  * x1[i] + beta2 *x2[i]+ b0[id[i]]
    predict[i]  ~ dlnorm(mu[i], sigma)C(lims[i,1],lims[i,2])
  }
  for (i in (N1+1):(N1+N2)){
    mu[i] <- beta0 + beta1  * x1[i] + beta2 *x2[i]+ b0[id[i]]
    predict[i]  ~ dlnorm(mu[i], sigma)C(lims[i,1],)
  }
  #priors
  sigma ~ dgamma(0.1, 0.1)
  tau_b0 <- 1/sigma_b0
  sigma_b0 ~ dunif(0,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
    b0[i] ~  dnorm(0,tau_b0)
  }
  for (i in 1:(N1+N2)){
    ppo[i] <- pow(2*3.141593,-0.5)*pow(sigma,-1)*exp(-0.5*pow((predict[i]-mu[i])/sigma, 2))
  }
}
write.model(model.function, "Scripts\\Bayes1_lognormal_Random.txt")
model.inits <- function(){list(sigma=2, beta0=1, beta1 = 1,beta2 = 1, b0 = c(rep(1,times = 20)), predict = Grub$lowerlim+1 )}
parameters = c("sigma", "beta2", "beta0", "beta1", "b0","ppo")


model.out <- bugs(model.data, model.inits, 
                  model.file = "Bayes1_lognormal_Random.txt",
                  parameters=parameters,
                  n.chains = 2,n.thin = 10, n.iter = 5000,  n.burnin = 1000, debug = T,
                  codaPkg=T,
                  working.directory = ".\\Scripts")

Lognormal_bayes <- read.bugs(model.out)



subset_pred <- grepl("ppo\\[", dimnames(Lognormal_bayes[[1]])[[2]])
mcmc_subset <- get_values(Lognormal_bayes,subset_pred)

#cpo which are now far of?
biggest <- abs(1/apply(as.matrix(mcmc_subset),2,mean))
plot(biggest)



subset_stuff <- grepl("predict\\[", dimnames(Lognormal_bayes_summary$statistics)[[1]])

Example123 <- Lognormal_bayes[[1]][123,subset_stuff]
Lognormal_bayes_summary <- summary(Lognormal_bayes)

predict_lognormal <- Lognormal_bayes_summary$statistics
#keine Ahnung was das jetzt aber bedeutet, wie weit sind die Werte von einer
#zufälligen Beoabchtung an der Stelle weg?
Example123[Example123>12] <- 12
#since it is censored, it does not look to normal 
hist(Grub$value - Example123, breaks = 100)


predict_cols <- grepl("predict\\[", dimnames(Lognormal_bayes_summary$statistics)[[1]])

Lognormal_bayes[1,predict_cols]



