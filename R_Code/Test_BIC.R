#Testfile BIC
source("helpfunctions.r")
library("runjags")
library("coda")
library("rjags")
set.seed(1341234234)
x <- runif(200,0,10)
y <- x + rnorm(200,0,5)

model.data <- list( y = y, N = length(y), x1 = x)
#DEFINE INTITIAL VALUES
model.inits <- list(list(sigma=2, beta0=1, beta1 = 1),
                    list(sigma=2, beta0=1, beta1 = 1)
)
#Monitored Variables
parameters <-c("beta0", "beta1", "sigma","Devi","D2evi")
#sigma is variance

model.function <- "model{
  for (i in 1:N){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- beta0 + beta1 *x1[i]
     D[i] <- - log(tau) + log(2*3.14159265358979) + pow(y[i]-mu[i],2)*tau
     D2[i] <- log(dnorm(y[i],mu[i],tau))
  }
  #priors
  Devi <- sum(D[])
  D2evi <- sum(D2[])
  tau <- 1/sigma
  sigma ~ dgamma(0.001, 0.001)
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
}"


runjags.options(method = "rjparallel")
#Set Up Model
#Generate MCMC SAMpls
Model_test <- run.jags(model = model.function,
                      monitor = parameters, data = model.data,
                      inits = model.inits, burnin = 2000,
                      sample = 5000, thin = 1, n.chains = 2)
plot(Model_test)
print(Model_test)
#The DIC Value for model comparison
dic_val <- extract.runjags(Model_test, "dic")
dic_val

mcmc_rep <- as.mcmc.list(Model_test)
subset_pred <- grepl("Devi", dimnames(mcmc_rep[[1]])[[2]])
mcmc_subset <- get_values(mcmc_rep,subset_pred)
mean(mcmc_subset)
subset_pred <- grepl("D2evi", dimnames(mcmc_rep[[1]])[[2]])
#works
mcmc_subset2 <- get_values(mcmc_rep,subset_pred)
-2 *mean(mcmc_subset2)
