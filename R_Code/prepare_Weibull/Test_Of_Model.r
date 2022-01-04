
#wenn nicht symmetrisch, kann es ja noch durch den random intercept erkläret werdne
#stimmt es zwar auch noc nicht ganz genau, soltle aber eingigermaßen passen


#plates
x2 <- rnorm(100, sd = 4)
x2 <- rep(x2,each = 100)
id <- as.factor(rep(1:100, each = 100))

#group
x3 <- rep(rnorm(20,sd = 2), each = 500)
group <- as.factor(rep(1:20, each = 500))
alkohol <- rnorm(10000, mean = 4, sd = 20)
y <- 3 + 3*alkohol + x2 + x3 + rnorm(100, sd = 10)

library(lme4)
summary(lm(y ~ id + group))

summary(lmer(y ~ group +alkohol+ (1|id)))
