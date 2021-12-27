#Bayes1_lognormal_random.r
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(coda)
library(R2OpenBUGS)

Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")
# Grub1 <- Grub
Grub <- Grub[Grub$group == 1,]
model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    id = (Grub$id)-10, Nsubj = length(unique(Grub$id)))

# MODEL SPECIFICATION 


model.function <- function(){
  for (i in 1:N){
    y[i] ~ dlnorm(mu[i], sigma)
    mu[i] <-  beta1 *x1[i] + b0[id[i]]
    predict[i]  ~ dlnorm(mu[i], sigma)
  }
  # forecast <- y[] - predict[]
  #priors
  sigma ~ dunif(0.1,100)
  beta1 ~ dnorm(0,0.000001)
  for ( i in 1:Nsubj){
    b0[i] ~ dnorm(0,0.0001)
  }
  
  
  sds <- pow(sd(b0[]),-1)
  for ( i in 1:Nsubj){
    # Distribution of future b0_i
    mu_gr[i] <- mean(mu[i:(i+6)])
    
    b0.rep[i] ~ dnorm(mu_gr[i],sds)  
  }
  
}
write.model(model.function, "Scripts\\Bayes1_lognormal_Random.txt")
model.inits <- function(){list(sigma=2, beta1 = 1, b0 = c(rep(1,times = 10)) )}
parameters = c("sigma", "beta1", "b0","predict","b0.rep")


model.out <- bugs(model.data, model.inits, 
                  model.file = "Bayes1_lognormal_Random.txt",
                  parameters=parameters,
                  n.chains = 2,n.thin = 10, n.iter = 5000,  n.burnin = 1000, debug = T,
                  codaPkg=T,
                  working.directory = ".\\Scripts")

Lognormal_bayes <- read.bugs(model.out)
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