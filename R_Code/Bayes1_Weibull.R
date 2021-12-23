#Weibull Modell 1 no random effects
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(coda)
library(R2OpenBUGS)
Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")

model.data <- list( y = Grub$value, N = length(Grub$value), x1 = Grub$grubsize,
                    x2 = Grub$group)



# MODEL SPECIFICATION 
#now the model specification is in line with the estimation in survival
model.function <- function(){
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
  k ~ dunif(0.1,100)
  beta0 ~ dnorm(0,0.000001)
  beta1 ~ dnorm(0,0.000001)
  beta2 ~ dnorm(0,0.000001)
}
write.model(model.function, "Scripts\\Bayes1_Weibull.txt")
model.inits <- function(){list(k=2, beta0=1, beta1 = 1,beta2 = 1 )}
parameters = c("scale", "beta2", "beta1", "beta0","predict")


model.out <- bugs(model.data, model.inits, 
                  model.file = "Bayes1_Weibull.txt",
                  parameters=parameters,
                  n.chains = 1, n.iter = 5000, n.thin = 4, n.burnin = 1000, debug = T,
                  codaPkg=T,
                  working.directory = ".\\Scripts")

Weibull_bayes <- read.bugs(model.out)
#k = 1.6
Weibull_summary <- summary(Weibull_bayes) #paramters are not like in the aft model bur rather as in the Ph (prop . hazard), therefore fliped (or it is the opposite direction what is more likely)
Weibull_summary$statistics





#parametric 
library(survival)
estimates <- survreg(Surv(value) ~ grubsize +  group, Grub, dist = "weibull")
estimates
1/0.62 # same in the weibull case if you flip the sign
Weibull_summary$statistics



#plot(survfit(Surv(value) ~ 1, data = Grub), 
plot(survfit(Surv(value) ~ group, data = Grub), 
     xlab = "Days", 
     ylab = "Overall survival probability",conf.int = F,col = c(1,2))
#For intervalls use Greenwalds formula