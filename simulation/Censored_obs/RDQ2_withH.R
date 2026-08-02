library(gamlss)
library(gamlss.cens)
library(VGAM)

RDQ2 <- function (mu.link="log", sigma.link="log", nu.link="log") {
    mstats <- checklink("mu.link", "Reparametrized.Dagum1", substitute(mu.link), c("inverse", "log", "identity"))# dummy
    dstats <- checklink("sigma.link", "Reparametrized.Dagum1", substitute(sigma.link), c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "Reparametrized.Dagum1", substitute(nu.link), c("inverse", "log", "identity"))

    structure(
          list(family = c("RDQ2", "Reparametrized.Dagum1"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), 
                nopar = 3,
                 type = "Continuous",
              mu.link = as.character(substitute(mu.link)), 
           sigma.link = as.character(substitute(sigma.link)), 
              nu.link = as.character(substitute(nu.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun,
           nu.linkfun = vstats$linkfun, 
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
           nu.linkinv = vstats$linkinv,
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                nu.dr = vstats$mu.eta,
                 dldd = function(y, mu, sigma, nu) {
                    tau = 0.5
                    (sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)/(sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)
                },
                 dldv = function(y, mu, sigma, nu) {
                    tau = 0.5
                    sigma*(-(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) - (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) - (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1))/(nu*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)
                },
                 dldm = function(y, mu, sigma, nu) {
                    tau = 0.5
                    sigma*(nu/mu)^sigma*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))/(mu*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)
                },
               d2ldd2 = function(y, mu, sigma, nu) {
                tau = 0.5
                (-sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)*log(y/nu)*log((nu/mu)^sigma + 1) - 2*sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)*log(nu/mu) + sigma*((y/nu)^sigma + 1)*(-sigma*(nu/mu)^(2*sigma)*((y/nu)^sigma + 1)*log(nu/mu)^2 - (nu/mu)^sigma*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)*log(nu/mu) + ((nu/mu)^sigma + 1)*(-sigma*(y/nu)^sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau)*log(y/nu)*log(nu/mu) + sigma*(y/nu)^sigma*(nu/mu)^sigma*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*(y/nu)^sigma*(nu/mu)^sigma*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log(nu/mu) - 2*sigma*(y/nu)^sigma*(nu/mu)^sigma*log(y/nu)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)^2*log((nu/mu)^sigma + 1) + sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*log(tau)*log(y/nu)*log(nu/y)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu)^2 + sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(tau)*log(nu/mu)*log(nu/y)*log((nu/mu)^sigma + 1) + sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(tau)*log(nu/mu)*log(nu/y) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)^2*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*((y/nu)^sigma*log(y/nu) + ((y/nu)^sigma + 1)*log(nu/y))*log(tau)*log(nu/mu) + (y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) + (y/nu)^sigma*((nu/mu)^sigma + 1)*log(y/nu)*log((nu/mu)^sigma + 1)^2 - (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1)^2 + (nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1)))*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*(-sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) + sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) + sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) - sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) - ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)*log((nu/mu)^sigma + 1))/(sigma^2*((y/nu)^sigma + 1)^2*((nu/mu)^sigma + 1)^2*log((nu/mu)^sigma + 1)^3)
            },
              d2ldddv = function(y, mu, sigma, nu) {
                tau = 0.5
                ((y/nu)^sigma*((nu/mu)^sigma + 1)*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)*log((nu/mu)^sigma + 1) - (nu/mu)^sigma*((y/nu)^sigma + 1)*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)*log((nu/mu)^sigma + 1) - 2*(nu/mu)^sigma*((y/nu)^sigma + 1)*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2) + ((y/nu)^sigma + 1)*(-sigma*(nu/mu)^(2*sigma)*((y/nu)^sigma + 1)*log(nu/mu) - sigma*(nu/mu)^sigma*((nu/mu)^sigma + 1)*log(tau)*log(nu/mu) + ((nu/mu)^sigma + 1)*(sigma*(y/nu)^sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau)*log(nu/mu) + sigma*(y/nu)^sigma*(nu/mu)^sigma*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) + sigma*(y/nu)^sigma*(nu/mu)^sigma*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu) - sigma*(y/nu)^sigma*(nu/mu)^sigma*log(y/nu)*log((nu/mu)^sigma + 1) + sigma*(y/nu)^sigma*(nu/mu)^sigma*log(nu/mu)*log((nu/mu)^sigma + 1) - sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) + sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(tau)*log(nu/y) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) - (y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) - (y/nu)^sigma*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2 - (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1)^2 + (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1)))*log((nu/mu)^sigma + 1))/(nu*((y/nu)^sigma + 1)^2*((nu/mu)^sigma + 1)^2*log((nu/mu)^sigma + 1)^3)
            },
              d2ldmdd = function(y, mu, sigma, nu) {
                tau = 0.5
                (nu/mu)^sigma*(-sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))*log(nu/mu)*log((nu/mu)^sigma + 1) - 2*sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))*log(nu/mu) + sigma*((nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu) + ((nu/mu)^sigma + 1)*((y/nu)^sigma*log(y/nu) + ((y/nu)^sigma + 1)*log(nu/y))*log(tau))*log((nu/mu)^sigma + 1) + (sigma*log(nu/mu) + 1)*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1))/(mu*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)^2*log((nu/mu)^sigma + 1)^3)
            },
               d2ldv2 = function(y, mu, sigma, nu) {
                tau = 0.5
                sigma*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(-(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) - (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) - (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) + sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*((y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) + (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1) - ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) + 2*sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*((y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) + (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1) - ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1)) + sigma*((y/nu)^sigma + 1)*(-(nu/mu)^(2*sigma)*((y/nu)^sigma + 1) - (nu/mu)^sigma*((nu/mu)^sigma + 1)*log(tau) + ((nu/mu)^sigma + 1)*((y/nu)^sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) - (y/nu)^sigma*(nu/mu)^sigma*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) - (y/nu)^sigma*(nu/mu)^sigma*(log(tau) - log((nu/mu)^sigma + 1)) + 2*(y/nu)^sigma*(nu/mu)^sigma*log((nu/mu)^sigma + 1) + (y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) - (y/nu)^sigma*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1) - (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log(tau) - (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1)))*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*((y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) + (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1) - ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1))/(nu^2*((y/nu)^sigma + 1)^2*((nu/mu)^sigma + 1)^2*log((nu/mu)^sigma + 1)^3)
            },
              d2ldmdv = function(y, mu, sigma, nu) {
                tau = 0.5
                sigma^2*(nu/mu)^sigma*(-(nu/mu)^sigma*((y/nu)^sigma + 1)*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) - 2*(nu/mu)^sigma*((y/nu)^sigma + 1)*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1)) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) + ((nu/mu)^sigma*((y/nu)^sigma + 1) + ((nu/mu)^sigma + 1)*log(tau))*log((nu/mu)^sigma + 1))/(nu*mu*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)^2*log((nu/mu)^sigma + 1)^3)
            },
            d2ldm2 = function(y, mu, sigma, nu) {
                tau = 0.5
                sigma*(nu/mu)^sigma*(2*sigma*(nu/mu)^sigma*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1)) + sigma*(nu/mu)^sigma*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1) - 1)*log((nu/mu)^sigma + 1) - ((nu/mu)^sigma + 1)*(sigma*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1)) + (sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1))/(mu^2*((nu/mu)^sigma + 1)^2*log((nu/mu)^sigma + 1)^3)
            },
          G.dev.incr  = function(y,mu,sigma,nu,...) -2*dRDQ2(y, mu, sigma, nu, log = TRUE), 
                rqres = expression(rqres(pfun="pRDQ2", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
           mu.initial = expression(mu <- rep(median(y),length(y))), 
        sigma.initial = expression(sigma <- rep(1, length(y))), 
           nu.initial = expression(nu <- rep(1, length(y))),
             mu.valid = function(mu) all(mu > 0), 
          sigma.valid = function(sigma)  all(sigma > 0),
             nu.valid = function(nu) all(nu > 0),
              y.valid = function(y)  all(y > 0),
                 mean = function(mu, sigma, nu) {
                    a = sigma
                    b = nu
                    p = -log(.5) / log((nu / mu)^sigma + 1)
                    return(ifelse(a > 1, p * b * beta(p+1/a, 1-1/a), Inf))
                  },
             variance = function(mu, sigma, nu) {
                a = sigma
                b = nu
                p = -log(.5) / log((nu / mu)^sigma + 1)
                return(ifelse(a > 2, b^2 * (beta(1-2/a, p+2/a) - beta(1-1/a, p+1/a)^2), Inf))
              }
          ),
                class = c("gamlss.family","family"))
}

dRDQ2 <- function(x, mu = 1, sigma = 1, nu = 1.5, log = FALSE) {
  tau <- 0.5
  mu2 <- -log(tau) / log((nu / mu)^sigma + 1)
  ddagum(x = x, scale = nu, shape1.a = sigma, shape2.p = mu2, log = log)
}

pRDQ2 <- function(q, mu = 1, sigma = 1, nu = 1.5, lower.tail = TRUE, log.p = FALSE) {
  tau <- 0.5
  mu2 <- -log(tau) / log((nu / mu)^sigma + 1)
  pdagum(q = q, scale = nu, shape1.a = sigma, shape2.p = mu2, lower.tail = lower.tail, log.p = log.p)
}

qRDQ2 <- function(p, mu = 1, sigma = 1, nu = 1.5, lower.tail = TRUE, log.p = FALSE) {
  tau <- 0.5
  mu2 <- -log(tau) / log((nu / mu)^sigma + 1)
  qdagum(p = p, scale = nu, shape1.a = sigma, shape2.p = mu2, lower.tail = lower.tail, log.p = log.p)
}

rRDQ2 <- function(n, mu = 1, sigma = 1, nu = 1.5) {
  tau <- 0.5
  mu2 <- -log(tau) / log((nu / mu)^sigma + 1)
  rdagum(n = n, scale = nu, shape1.a = sigma, shape2.p = mu2)
}

# Log Verosimilitud

library(pracma)
LogLS <- function(theta, DM, data) {
  mu_ <- theta[1:3]; sigma_ <- theta[4:5]; nu_ <- theta[6]
  
  X1 <- DM[, 1:3]; X2 <- DM[, 4:5]; X3 <- DM[, 6]

  mu    <- as.vector(exp(X1 %*% mu_))       # Log-link
  sigma <- as.vector(exp(X2 %*% sigma_))    # Log-link
  nu    <- as.vector(exp(X3 * nu_))         # Log-link

  T <- data[, 1]
  I <- data[, 2]
  indx_C <- which(I == 0)

  T_C <- T[indx_C]
  mu_C <- mu[indx_C]; sigma_C <- sigma[indx_C]; nu_C <- nu[indx_C]
  
  T_F <- T[-indx_C]
  mu_F <- mu[-indx_C]; sigma_F <- sigma[-indx_C]; nu_F <- nu[-indx_C]

  return(sum(dRDQ2(T_F, mu = mu_F, sigma = sigma_F, nu = nu_F, log = TRUE)) +
         sum(pRDQ2(T_C, mu = mu_C, sigma = sigma_C, nu = nu_C, lower.tail = FALSE, log.p = TRUE)))
}

#######################################
############# Monte Carlo #############
#######################################

library(parallel)
source("./find_delta.R")

MCllikelihood_estimation <- function(n = 1000, mu_ = c(1, .5, .2), sigma_ = c(.5, .4), nu_ = c(-.2), censorship = 0.1) {

  z1 <- as.numeric(scale(runif(n)))
  z2 <- as.numeric(scale(rnorm(n)))
  
  X1 <- model.matrix(~z1+z2)                   #matriz de diseño de mu
  X2 <- model.matrix(~z1)                      #matriz de diseño de sigma
  X3 <- matrix(1, nrow = n)                    #matriz de diseño de nu
  colnames(X3) <- "(Intercept)"
  
  mu.true <- as.vector(exp(X1 %*% mu_))        # log-link
  sigma.true <- as.vector(exp(X2 %*% sigma_))  # log-link
  nu.true <- as.vector(exp(X3 %*% nu_))        # log-link-desplazado (nu > 1)
  y <- rRDQ2(n, mu.true, sigma.true, nu.true)  # simulando valores
  lmbds <- apply(cbind(mu.true, sigma.true, nu.true), 1, find_delta,
                 model = dRDQ2, censorship = censorship, search_interval = c(0, 100))
  C <- rexp(n = n, lmbds)
  
  # GAMLSS
  aux <- gamlss(Surv(time = pmin(y, C), event = y < C) ~ X1[, -1, drop=FALSE],
                sigma.fo = ~X2[, -1, drop = FALSE],
                family = cens(RDQ2), method = RS(1000),
                control = gamlss.control(trace = FALSE))
  capture.output(res.gamlss <- summary(aux)[, 1:2], file = nullfile())

  res <- cbind(c(mu_, sigma_, nu_),                              # REAL
               res.gamlss)                                       # GAMLSS

  H <- hessian(LogLS, x0 = res[, 2], DM = cbind(X1, X2, X3),
               data = Surv(time = pmin(y, C), event = y < C))

  VCOV <- sqrt(diag(solve(-H)))

  res <- cbind(res, VCOV)

  colnames(res) <- c("true",
                     "est gamlss", "se gamlss", "se hessian")
  rownames(res) <- c(paste("beta1", 1:ncol(X1), sep = ""),
                     paste("beta2", 1:ncol(X2), sep = ""),
                     paste("beta3", 1:ncol(X3), sep = ""))
  return(list(Results = res, LL = c(gamlss = logLik(aux)), Converged = aux$converged))
}

# to use in mclapply
RNGkind("L'Ecuyer-CMRG")

set.seed(123)
seeds <- list(.Random.seed)
for (i in 2:4000) {
  seeds[[i]] <- nextRNGStream(seeds[[i - 1]])
}

prob <- function(i, theta, seeds, casos.n, censorship) {
  if (i <= 1000) {
    n <- casos.n[1]
  } else if (i <= 2000) {
    n <- casos.n[2]
  } else if (i <= 3000) {
    n <- casos.n[3]
  } else {
    n <- casos.n[4]
  }

  .Random.seed <<- seeds[[i]]
  temp <- TRUE
  errors <- -1
  while (temp) {
    testing <- try(MCllikelihood_estimation(n = n,
                                            mu_ = theta[1:3],
                                            sigma_ = theta[4:5],
                                            nu_ = theta[6],
                                            censorship = censorship),
                   silent = TRUE)
		temp <- grepl("Error", testing)[1]
    errors <- errors + 1
  }
  testing[["Errors"]] <- errors
  cat("%:", round(i / 4000 * 100, 3), " (", i, "/", 4000, ")  n:", n, "\n", sep = "")
  return(testing)
}

# sim

casos.n <- c(50, 100, 200, 500)
casos.par <- matrix(c(1, .5, .2, .5, .4, -.2,
                      -1, -.5, -.2, .25, .7, .3),
                    nrow = 2, byrow = TRUE)

# mclapply solo funciona en Linux y MacOS,
# no en windows. Si corre esto en windows hace un lapply solamente.
# Forma mas facil de aplicar calculo paralelo.
for (k in 3) {

  Cens <- c(.1, .25, .5)[k]

  Cens_level <- c("10", "25", "50")[k]

  system.time(testing.p1 <- mclapply(1:4000, prob,
                                     seeds = seeds,
                                     casos.n = casos.n,
                                     theta = casos.par[1, ],
                                     censorship = Cens)) # 8 min approx
  system.time(testing.p2 <- mclapply(1:4000, prob,
                                     seeds = seeds,
                                     casos.n = casos.n,
                                     theta = casos.par[2, ],
                                     censorship = Cens)) # 7 min approx
  
  setwd(paste("./RDQ2/RDQ2_", Cens_level, sep = ""))
  
  reps <- 1000
  for (i in 1:length(casos.n)) {
    
    ############ set 1
  
    temp_p1 <- testing.p1[(1 + reps * (i - 1)):(reps * i)]
  
    # LL
    LL_p1 <- sapply(temp_p1, "[[", 2)
  	name_ll_p1 <- paste("LL_p1_", casos.n[i], '.csv', sep = '')
    write.csv(LL_p1, file = name_ll_p1, row.names = FALSE)
  
    # pars
    pars_p1 <- t(sapply(temp_p1, "[[", 1))
  	name_pars_p1 <- paste("pars_p1_", casos.n[i], '.csv', sep = '')
    write.csv(pars_p1, file = name_pars_p1, row.names = FALSE)
    
    # Converged
    conv_p1 <- sapply(temp_p1, "[[", 3)
  	name_conv_p1 <- paste("noconv_p1_", casos.n[i], '.csv', sep = '')
    write.csv(sum(!conv_p1), file = name_conv_p1, row.names = FALSE)
    
    # Errors
    errors_p1 <- t(sapply(temp_p1, "[[", 4))
  	name_errors_p1 <- paste("errors_p1_", casos.n[i], '.csv', sep = '')
    write.csv(sum(errors_p1), file = name_errors_p1, row.names = FALSE)
  
    ############ set 2
  
    temp_p2 <- testing.p2[(1 + reps * (i - 1)):(reps * i)]
  
    # LL
    LL_p2 <- sapply(temp_p2, "[[", 2)
  	name_ll_p2 <- paste("LL_p2_", casos.n[i], '.csv', sep = '')
    write.csv(LL_p2, file = name_ll_p2, row.names = FALSE)
  
    # pars
    pars_p2 <- t(sapply(temp_p2, "[[", 1))
  	name_pars_p2 <- paste("pars_p2_", casos.n[i], '.csv', sep = '')
    write.csv(pars_p2, file = name_pars_p2, row.names = FALSE)
    
    # Converged
    conv_p2 <- sapply(temp_p2, "[[", 3)
  	name_conv_p2 <- paste("noconv_p2_", casos.n[i], '.csv', sep = '')
    write.csv(sum(!conv_p2), file = name_conv_p2, row.names = FALSE)
  
    # Errors
    errors_p2 <- t(sapply(temp_p2, "[[", 4))
  	name_errors_p2 <- paste("errors_p2_", casos.n[i], '.csv', sep = '')
    write.csv(sum(errors_p2), file = name_errors_p2, row.names = FALSE)
  }

}
