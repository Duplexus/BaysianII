#Help FUnctios
Grub <- read.csv("../data/Grubs_Easy_normalized_size.csv")


# #parametric
# library(survival)
# lognormal_results1 <- survreg(Surv(value) ~ grubsize +  group, Grub, dist = "lognormal")
# # scale is sigma seems to be 1/the one from the baysian but beside this similar
# #results between bayes and model are really similar
# summary(lognormal_results1)
# 
# ## Model is runned with predictions again
# ####PPC####
# #wie sieht das gesch?tzte Modell aus:
# subset_pred <- grepl("predict\\[", dimnames(lognormal1[[1]])[[2]])
# #extracts all values and packages them into one list
# #... one can include which rows are wanted
# #object is object after coda
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

# data_conv <- (get_values(lognormal1,subset_pred))
# #mean similar?
# vgl_fct(data_conv,fct = mean, mean(Grub$value))
# #sd similar
# vgl_fct(data_conv,fct = sd, sd(Grub$value))

#wie viel Prozent der Residuuen sind ?ber zwei Standardeinheiten weg?
#takes a vector
outer_5 <- function(x,sds=1.96){
  l <- length(x)
  a1 <- x >= sds*sd(x) + mean(x)
  a2 <- x <= - sds*sd(x) + mean(x)
  return(list(left = sum(a2)/l, right = sum(a1)/l, comb = sum(a2)/l+ sum(a1)/l))
  
}


#### lognorm_lognormal_random ####
# lognorm_dic(lognorm_rand_cens_mcmc_rep)

lognorm_dic <- function(lognorm_object){
  lognorm_object = as.mcmc.list(lognorm_object)
  lognorm_object2 = as.mcmc.list(lognorm_rand_cens_sens4)
  
  
  subset_pred <- grepl("Deviance", dimnames(lognorm_object[[1]])[[2]])
  lognorm_subset <- get_values(lognorm_object,subset_pred)
  md <- mean(lognorm_subset)

  a <- vector()
  subset_pred <- grepl("beta0", dimnames(lognorm_object[[1]])[[2]])
  a["beta0"] <- mean(get_values(lognorm_object,subset_pred))

  subset_pred <- grepl("beta2", dimnames(lognorm_object[[1]])[[2]])
  a["beta2"] <- mean(get_values(lognorm_object,subset_pred))

  subset_pred_b1 <- grepl("beta1", dimnames(lognorm_object[[1]])[[2]])
  if (any(subset_pred_b1)){
  a["beta1"] <- mean(get_values(lognorm_object,subset_pred_b1))
  }
  subset_pred <- grepl("^tau$", dimnames(lognorm_object[[1]])[[2]])
  a["tau"] <- mean(get_values(lognorm_object,subset_pred))
  a1 <- a
  # cat("order beta0 beta2 beta1 tau\n")
  # cat("means:: ",a,"\n\n")

  subset_pred <- grepl("b0\\[", dimnames(lognorm_object[[1]])[[2]])
  b0_subset1 <- apply(get_values(lognorm_object,subset_pred),2,mean)
  b0_subset <- b0_subset1[Grub$id]

  subset_pred <- grepl("y\\[", dimnames(lognorm_object[[1]])[[2]])
  lognorm_subset_y <- get_values(lognorm_object,subset_pred)

    Grub$sim_value <- apply(lognorm_subset_y,2,mean)
  if (any(subset_pred_b1)){
    pd <- md - (-2 *sum(log(dlnorm(Grub$sim_value,a["beta0"]+a["beta1"]*Grub$grubsize
                                   + a["beta2"]*Grub$group + b0_subset,sqrt(1/a["tau"])))))
  }else{
    pd <- md - (-2 *sum(log(dlnorm(Grub$sim_value,a["beta0"]
                                   + a["beta2"]*Grub$group + b0_subset,sqrt(1/a["tau"])))))
  }
  cat("pd: \t",pd,"\nDeviance:", md,"\nDIC: \t",pd+md,"\n")
  return(c(pd,md,pd+md))
}
