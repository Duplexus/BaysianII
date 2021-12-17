setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
 
library(dplyr)
library(ggplot2)
Grub <- read.csv("..\\data\\Grubs_Nematodes.csv")
NAs <- is.na(Grub$UPPERLIM)
Grub$UPPERLIM[NAs]<- 12
Grub$value <- (Grub$UPPERLIM + Grub$LOWERLIM) /2
Grub$UPPERLIM[NAs]<- NA
rm(NAs)
colnames(Grub) <- c("group", "id", "grubsize", "lowerlim", "upperlim", "value") 
Grub$group <- as.factor(Grub$group)
#x ~ dweib(v, lambda)

Grub %>% group_by(group) %>% summarise(across(where(is.numeric),
                                       function(x)sum(is.na(x))))
#mean
Grub %>% group_by(group) %>% summarise(across(where(is.numeric),function(x) mean(x, na.rm  = T)))
#group 1 seems to survive longer
#median
Grub %>% group_by(group) %>% summarise(across(where(is.numeric),function(x) median(x, na.rm  = T)))

linmod <- lm(value ~ group + grubsize, Grub)
summary(linmod)
#glinmod <- glm(value ~ group + grubsize, Grub,family =poisson())

glm.fit()

