#Masterfile:
source("helpfunctions.r")
library("runjags")
library("coda")
library("rjags")

####Without censored, without random effects####
# Lognormal----
source("Bayes1_lognormal_jags.R")
plot(lognorm)
print(lognorm)
#The DIC Value for model comparison
dic_val <- extract.runjags(lognorm, "dic")
dic_val
extract.runjags(ModelLogN, "stochastic")

#coda
mcmc_lognorm <- as.mcmc.list(lognorm)
#convergence
gelman.diag(mcmc_lognorm, confidence = 0.95)
gelman.plot(mcmc_lognorm, confidence = 0.95)
geweke.diag(mcmc_lognorm)
geweke.plot(mcmc_lognorm)
plot(mcmc_lognorm)

