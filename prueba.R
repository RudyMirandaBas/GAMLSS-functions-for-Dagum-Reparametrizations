library(dplyr)
library(gamlss)

source("./RD.R")
source("./RDQ1.R")
source("./RDQ2.R")

z1 <- scale(rnorm(100))
z2 <- scale(runif(100))

mu_coef <- c(2, -1, 1)
sigma_coef <- c(-1, .5)
nu_coef <- 0.5

mu <- as.vector(exp(model.matrix(~ z1 + z2) %*% mu_coef))
#sigma <- as.vector(exp(model.matrix(~ z1) %*% sigma_coef))
sigma <- rep(exp(-1), 100)
nu <- rep(exp(nu_coef) + 1, 100)

y <- rRD(100, mu, sigma, nu)

# checking 
dagum <- gamlss(y ~ z1 + z2, sigma.fo = ~ 1,
								family = GB2, method = RS(100),
								tau.fix = TRUE, tau.start = 1)
dagum = gamlss(y ~ z1 + z2, family = WEI)
rd <- gamlss(y ~ z1 + z2, sigma.fo = ~ 1,
						 family = WEI3, method = RS(100))
TAU <- .5
rdq1 <- gamlss(y ~ z1 + z2, sigma.fo = ~ z1,
							 family = RDQ1, method = RS(100))
rdq2 <- gamlss(y ~ z1 + z2, family = RDQ2, method = RS(100))

list(dagum = logLik(dagum), rd = logLik(rd),
		 rdq1 = logLik(rdq1), rdq2 = logLik(rdq2))








d1 <- dWEI3(y, mu = fitted(rd, "mu"),
          sigma = fitted(rd, "sigma"),
          log = TRUE)

d2 <- dWEI(y, mu = fitted(dagum, "mu"),
           sigma = fitted(dagum, "sigma"),
           log = TRUE)

max(abs(d1 - d2))


source("https://git.io/JRi5k")
source("./../DATA/RBS.R")

library(RelDist)
