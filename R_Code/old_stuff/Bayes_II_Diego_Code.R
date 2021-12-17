# Bayesian II - Project Part I

#############################################
# Import
#############################################
library(coda)
library(ggplot2)
library(R2OpenBUGS)

setwd("C:/Users/Diego/OneDrive - KU Leuven/KUL/3_sem/Bayes II")
data <- read.csv("Grubs_Nematodes.csv")

# modidy data
data1<-data;
data1[4][is.na(data1[5])] # all lower lims are 12
data1[5][is.na(data1[5])] <- 12;

#############################################
# 1: Median death time
#############################################
# data1['MEAN'] = rowMeans(data1[c('LOWERLIM','UPPERLIM')])
data1['MIDPOINT'] = apply(data1[c('LOWERLIM','UPPERLIM')], 1, median)
data1['Center_size'] = data1['GRUBSIZE'] - apply(data1['GRUBSIZE'], 2, mean)
data1['std_size'] = scale(data1['GRUBSIZE'])[1:140]

group1 = data1$MIDPOINT[data1['GROUP']==1]
group2 = data1$MIDPOINT[data1['GROUP']==2]
median(group1) # =5
median(group2) # =3
#graph
data1$GROUP = as.factor(data1$GROUP)
g1<-ggplot(data1, aes(x=MIDPOINT, color=GROUP))
  g1 + 
  geom_histogram(binwidth = 1, fill='white', position = 'identity', alpha=0.4) +
  geom_vline(aes(xintercept=c(3),color='2')) +
  geom_vline(aes(xintercept=c(5),color='1')) 
  
#############################################
# 2: Grubs' size covariate
#############################################
# effect of $GRUBSIZE on $MEAN: mean ~ beta*size
g2<-ggplot(data1, aes(x=GRUBSIZE, y=MIDPOINT, color=GROUP))
g2 + geom_point()
  
# DEFINE INITIAL VALUES
data1$GROUP = as.integer(data1$GROUP)
data1$GROUP = as.numeric(data1$GROUP==1)

model.inits <- list(beta0=rnorm(1), beta_G=rnorm(1), beta=rnorm(1), sigma=rnorm(1))
model.data <- list(
  # x = data1$GRUBSIZE, g = data1$GROUP, y = data1$MIDPOINT, N = nrow(data1)
  x = data1$Center_size, g = data1$GROUP, y = data1$MIDPOINT, N = nrow(data1)
  # x = data1$std_size, g = data1$GROUP, y = data1$MIDPOINT, N = nrow(data1)
)

# MODEL SPECIFICATION 
survival1 <- function(){
  # Specification data model
  for (i in 1:N)
  {
    y[i] ~ dnorm(y.hat[i],tau)
    y.hat[i]<- (beta0 + beta_G*g[i] + beta*x[i])
  }
  # Prior specification
  beta ~ dunif(-100,100) 
  beta0 ~ dunif(-100,100)
  beta_G ~ dunif(-100,100)
  # beta ~ dunif(-10,10) # not working, takes too long?
  # beta0 ~ dunif(-10,10)
  tau <- pow(sigma, -2)
  sigma ~ dunif(0,2)
  # sigma ~ dnorm(0.5,0.25)
}
write.model(survival1, "survival1.txt")

# these are the parameters to save
parameters = c("beta0", "beta_G","beta", "sigma")

# SET UP AND RUN MODEL
# specify model, data, number of parallel chains
model.out <- bugs(model.data, model.inits, 
                  model.file = "survival1.txt",
                  parameters=parameters,
                  n.chains = 4, n.iter = 20000,  n.burnin = 10000, debug = T,
                  codaPkg=T, working.directory = getwd())

# CONVERGENCE:
#   - using x = grubsize, there is persistency in the traceplots beta0 and beta
#   - try centering/stdizing. Using x = center_size Convergence improves!

# posterior summary statistics
out <- read.bugs(model.out)

HPDinterval(as.mcmc(as.matrix(out)))

plot(out)
  
# Posterior distributions
# densplot(out) 

# Posterior summary statistics
summary(out)
#             2.5%     25%     50%      75%    97.5%
# beta      -2.6110  -1.816  -1.401  -0.9842  -0.1813
# beta0      2.9540   3.262   3.422   3.5830   3.8920
# beta_G     0.9213   1.356   1.583   1.8160   2.2520
# deviance 736.7000 738.600 740.200 742.4000 748.7000
# sigma      1.9530   1.982   1.991   1.9960   2.0000

#############################################
# 3: Random effect
#############################################
# effect of $UREPID on last model
g2<-ggplot(data1, aes(x=GRUBSIZE, y=MIDPOINT, color=GROUP))
g2 + geom_point()

# DEFINE INITIAL VALUES
model.inits <- list(beta=rnorm(1), betea_G=rnorm(1), beta0=rnorm(20), sigma=rnorm(1), 
                    sigmab=rnorm(20))
model.data <- list(
  x = data1$Center_size, g = data1$GROUP, y = data1$MIDPOINT, N = nrow(data1), id = data1$UREPID
  # x = data1$std_size, g = data1$GROUP, y = data1$MIDPOINT, N = nrow(data1), id = data1$UREPID
)

# MODEL SPECIFICATION
survival1mix <- function(){
  # Specification data model
  for (i in 1:N)
  {
    y[i] ~ dnorm(y.hat[i],tau)
    y.hat[i]<- (beta0[id[i]] + beta_G*g[i]  + beta*x[i])
  }
  
  for (j in 1:20)
  {
    beta0[j] ~ dnorm(0, taub0)
  }
  
  # Prior specification
  beta ~ dunif(-100,100)
  beta_G ~ dunif(-100,100)
  tau <- pow(sigma, -2)
  sigma ~ dunif(0,2)
  taub0 <- pow(sigmab, -2)
  sigmab ~ dunif(0,100)
}

write.model(survival1mix, "survival1mix.txt")

# these are the parameters to save
parameters = c("beta0","beta", "beta_G", "sigma","sigmab")

# SET UP AND RUN MODEL
# specify model, data, number of parallel chains
model.out <- bugs(model.data, model.inits, 
                  model.file = "survival1mix.txt",
                  parameters=parameters,
                  n.chains = 5, n.iter = 20000,  n.burnin = 10000, debug = T,
                  codaPkg=T,working.directory = getwd())

# CONVERGENCE:
#   - using x = center_size there is some autocorr for sigmab and sigma
#   - using x = std_size, same problem
#   - turning the over.relax option on, same problem.  over.relax = TRUE
#   - try two level hierarchical model?

# posterior summary statistics
out <- read.bugs(model.out)

HPDinterval(as.mcmc(as.matrix(out)))

plot(out) # see in openbugs

# Posterior distributions
densplot(out) # see in openbugs

# Posterior summary statistics
summary(out)
#              2.5%      25%      50%      75%    97.5%



#############################################
# 4: Outliers or Influential observations
#############################################


#############################################
# 5: Distribution of the Random effect
#############################################


#############################################
# 6: Interval censored?
#############################################
