#First Bayes test:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

getwd()
library(coda)
library(R2OpenBUGS)
Grub <- read.csv("..\\data\\Grubs_Nematodes.csv")
NAs <- is.na(Grub$UPPERLIM)
Grub$UPPERLIM[NAs]<- 12
Grub$value <- (Grub$UPPERLIM + Grub$LOWERLIM) /2
Grub$UPPERLIM[NAs]<- NA
rm(NAs)
colnames(Grub) <- c("group", "id", "grubsize", "lowerlim", "upperlim", "value") 
#wenn ich das als Faktor definere rafft das WIbugs nicht mehr
#Grub$group <- as.factor(Grub$group)
Grub$group <- Grub$group-1

#https://www.youtube.com/watch?v=_44_RXTWpRw
#lognormal with value
N <- nrow(Grub)
y <- log(Grub$value)
x1 <- Grub$grubsize
x2 <- Grub$group
data <- list ("N", "y", "x1", "x2")
inits <- function(){list(beta0 = 0, beta1 = 0, beta2 = 0, sigma = 3)
}
getwd()
# setwd(".\\..\\R_code")
# Grub_sim <- bugs(data, inits, model.file = "Bayes1.txt",codaPkg=T,
#                     parameters = c("beta0", "beta1","beta2","sigma"),
#                     n.chains = 3,n.burnin = 0, n.iter = 10000, n.thin = 10)#, debug=TRUE)
# # setwd("..")
# typeof(Grub_sim)
# Grub_sim_coda <- read.bugs(Grub_sim)
# typeof(Grub_sim_coda)
# pdf("Test1.pdf")
# geweke.plot(Grub_sim_coda, frac1=0.1, frac2=0.5)
# dev.off()



#how many obs per id 
Nsubj <- as.vector(table(Grub$id))
CumSum <- cumsum(Nsubj) - 7
N <- length(CumSum)
y <- log(Grub$value)
x1 <- Grub$grubsize
x2 <- Grub$group
z1 <- Grub$id
data2 <- list ("N","Nsubj","CumSum", "y", "x1", "x2", "z1")
subv <- rep(0.1,times = 20)
inits2 <- function(){list(beta0.sub = 0, beta1.sub = 0,
                          beta2.sub = 0,sigma.sub = 0.001, 
                          b0.grp = subv, sigma_b0.grp=subv)
}
bugs.data(data2, dir=getwd(),digits=5,data.file="data.txt")	
#beta0 = Intercept , beta1 = Grubsize, beta2 = Grub Group, b0 random intercept 
#per individuum
setwd(".\\R_code")
Grub_sim_random <- bugs(data2, inits2, model.file = "lognormal_try1.txt",
                        codaPkg=T,parameters = c("beta0.sub", "beta1.sub",
                                                 "beta2.sub","sigma.sub",
                                                 "b0.grp","sigma_b0.grp"),
                        n.chains = 1,n.burnin = 5000, n.iter = 10000, n.thin = 1, debug=T)
setwd("..")
typeof(Grub_sim_random)
Grub_sim_coda <- read.bugs(Grub_sim_random)
typeof(Grub_sim_coda)
# pdf("Test1.pdf")
# geweke.plot(Grub_sim_coda, frac1=0.1, frac2=0.5)
# dev.off()
summary(Grub_sim_coda)

#Try a LMM Model to counter it
library(lme4)
lmer()
colnames(Grub)
Grub$id <- as.factor(Grub$id)
lmm_model <- lmer(I(log(value)) ~ grubsize + group + (1|id),Grub)
summary(lmm_model)
glmm_model <- glmer(value ~ grubsize + group + (1|id),Grub,family = poisson(link = "log"))
summary(glmm_model)

library(corrplot)
library
Korrelationsmatrix <- crosscorr(Grub_sim_coda)
#relativ geringe Korrelation zwischen den Schätzern
corrplot(Korrelationsmatrix[21:25,21:25], method = 'color', order = 'alphabet')
autocorr.plot(Grub_sim_coda)
batchSE(Grub_sim_coda)
crosscorr.plot(Grub_sim_coda)
#cumuplot(Grub_sim_coda)
effectiveSize(Grub_sim_coda)
gelman.plot(Grub_sim_coda)
plot.mcmc(Grub_sim_coda)
traceplot(Grub_sim_coda,col = "red")


survreg()