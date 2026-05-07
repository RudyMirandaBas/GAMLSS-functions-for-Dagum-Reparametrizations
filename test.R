library(gamlss)

source("./RD_fixed.R")
source("./RDQ1.R"); TAU <- .5
source("./RDQ2.R"); TAU <- .5

sample <- rRD(1000, mu = 1, sigma = 2, nu = 2); mean(sample)

rd <- gamlss(sample ~ 1, family = RD, method = RS(100))
rdq1 <- gamlss(sample ~ 1, family = RDQ1, method = RS(100))
rdq2 <- gamlss(sample ~ 1, family = RDQ2, method = RS(100))
d <- gamlss(sample ~ 1, family = GB2, method = RS(1000), tau.fix = TRUE, tau.start = 1)

rd
rdq1
d

source("./RDQ1_test.R")
sample <- rRDQ1(n = 1000, mu = 1, sigma = 1, nu = 1, tau = .25)
quantile(sample, .25)

gamlss(sample ~ 1, family = RDQ1, tau.fix = TRUE, tau.start = .25)
