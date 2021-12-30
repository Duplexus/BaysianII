#Help FUnctios
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
Grub <- read.csv("..\\data\\Grubs_Easy_normalized_size.csv")


# #parametric
# library(survival)
# lognormal_results1 <- survreg(Surv(value) ~ grubsize +  group, Grub, dist = "lognormal")
# # scale is sigma seems to be 1/the one from the baysian but beside this similar
# #results between bayes and model are really similar
# summary(lognormal_results1)
# 
# ## Model is runned with predictions again
# ####PPC####
# #wie sieht das geschätzte Modell aus:
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

# 
# data_conv <- (get_values(lognormal1,subset_pred))
# #mean similar?
# vgl_fct(data_conv,fct = mean, mean(Grub$value))
# #sd similar
# vgl_fct(data_conv,fct = sd, sd(Grub$value))

#wie viel Prozent der Residuuen sind über zwei Standardeinheiten weg?
#takes a vector
outer_5 <- function(x,sds=1.96){
  l <- length(x)
  a1 <- x >= sds*sd(x) + mean(x)
  a2 <- x <= - sds*sd(x) + mean(x)
  return(list(left = sum(a2)/l, right = sum(a1)/l, comb = sum(a2)/l+ sum(a1)/l))
  
}
# #just log because we live in the log world
# data <- log(data_conv[200,]) - log(Grub$value)
# outer_5(data)$comb
# hist(data)