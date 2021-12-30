

x <- runif(10,5,10)
x2 <- rep(x,each = 10)
x3 <- rep(c(2,5), each = 50)
group <- as.factor(rep(1:2, each = 50))
id <- as.factor(rep(1:10, each = 10))
y <- 3 + x2 + x3 + rnorm(100)


summary(lm(y ~ id + group))
