

fx <- function(x){
  if (x >= 0 & x <= 1){
    return(0.6)
  }else if( x >= 2 & x <= 3){
    return(0.4)
  }else{
    return(0)
  }
}
#vectorized version
posterior1 <- function(x){
  return(sapply(x,fx))
}
x <- seq(-2,5, by = 0.01)
plot(x, posterior1(x), type = "l")





acceptance <- function(new, old, posterior){
  min(c(posterior(new) / posterior(old),1))
}
acceptance(2,1,posterior)

2:5
set.seed(4124131)
metropolis <- function(leng, star, posterior, sd){
  values <- c(star,c(1:leng))
  for (i in 2:length(values)){
    old <- values[i-1]
    #using the proposal
    #normal prior
    new <- values[i-1] + rnorm(1,mean = 0, sd = sd)
    #uniform prior
    #new <- values[i-1] + runif(1,-sd,sd)
    acceptan <- acceptance(new, old, posterior)
    if (acceptan > runif(1)){
      values[i] <- new
    }else{
      values[i] <- old
    }
  }
  return(values[2:length(values)])
}


results <- metropolis(20000, 0, posterior1, 0.25)
hist(results)
####nächste Sache die Wahrscheinlichkeiten stecken um####
#rauszubekommen, warum das wirklich so ist


mean(dnorm(1, mean = seq(0,1, by = 0.001) ,sd = 0.2831))

pnorm(-0.5, sd = 0.25)*2

mean(dnorm(0, mean = seq(0,1, by = 0.001) ,sd = 0.25))

z <- seq(0,1, by = 0.01)

p_steck_unten <- function(z){
  steck_unten_approx <- function(z){
    return(sum(dnorm(z, mean = seq(-10,0, by = 0.001) ,sd = 0.25)))
  }

  return(steck_unten_approx(z)/ (steck_unten_approx(0)*2))
}

# #die Wahrscheinlichkeit auf Untenstecken für 0 (ist ja ungefähr
# # 0.5 da über 4 sd sehr unwahrscheinlich) normalisieren
p_steck_oben <- function(z){
  steck_oben_approx <- function(z){
    return(sum(dnorm(z, mean = seq(1,11, by = 0.001) ,sd = 0.25)))
  }
  return(steck_oben_approx(z)/ (steck_oben_approx(1)*2))
}

# p_steck_oben <- function(z){
#   return(p_steck_unten(1 - z))
# }

p_steck_unten(0.1)
p_steck_oben(0.9)


#normalisieren auf 0.5 
plot(z,sapply(z, p_steck_unten), type = "l")
plot(z,sapply(z, p_steck_oben), type = "l")

#logisches Ergebnis, da ist quasi die Frage, was ist die durchschnitt
#liche dichte der einzelnen Verteilungsfunktionen. Da wir auf 0,1 arbeiten
# mit einer geringen sd. Ist für den mittelsten Wert die durcschnittliche Dichte
# ungefähr 0.95, da symetrische Wahrscheinlichkeit und das dann die FLäche
p_getroffen <- function(z){
  return(mean(dnorm(z, mean = seq(0,1, by = 0.0001) ,sd = 0.25)))
}
plot(z,sapply(z, p_getroffen) + sapply(z, p_steck_oben) + sapply(z, p_steck_unten) , type = "l")

a <-  sapply(z, p_getroffen) + sapply(z, p_steck_oben) + sapply(z, p_steck_unten)
a





