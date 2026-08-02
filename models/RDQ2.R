library(gamlss)
library(VGAM)

RDQ2 <- function (mu.link="log", sigma.link="log", nu.link="log", tau.link="identity") {
    mstats <- checklink("mu.link", "Reparametrized.Dagum2", substitute(mu.link), c("inverse", "log", "identity"))
    dstats <- checklink("sigma.link", "Reparametrized.Dagum2", substitute(sigma.link), c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "Reparametrized.Dagum2", substitute(nu.link), c("inverse", "log", "identity"))
    tstats <- checklink("tau.link", "Reparametrized.Dagum2", substitute(tau.link), c("inverse", "log", "identity"))

    structure(
          list(family = c("RDQ2", "Reparametrized.Dagum2"),
           parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
                nopar = 4,
                 type = "Continuous",
              mu.link = as.character(substitute(mu.link)), 
           sigma.link = as.character(substitute(sigma.link)), 
              nu.link = as.character(substitute(nu.link)), 
             tau.link = as.character(substitute(tau.link)), 
           mu.linkfun = mstats$linkfun, 
        sigma.linkfun = dstats$linkfun,
           nu.linkfun = vstats$linkfun, 
          tau.linkfun = tstats$linkfun, 
           mu.linkinv = mstats$linkinv, 
        sigma.linkinv = dstats$linkinv,
           nu.linkinv = vstats$linkinv,
          tau.linkinv = tstats$linkinv,
                mu.dr = mstats$mu.eta, 
             sigma.dr = dstats$mu.eta, 
                nu.dr = vstats$mu.eta,
               tau.dr = tstats$mu.eta,
                 dldd = function(y, mu, sigma, nu, tau) {
                    (sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)/(sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)
                },
                 dldv = function(y, mu, sigma, nu, tau) {
                    sigma*(-(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) - (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) - (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1))/(nu*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)
                },
                 dldm = function(y, mu, sigma, nu, tau) {
                    sigma*(nu/mu)^sigma*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))/(mu*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)
                },
               d2ldd2 = function(y, mu, sigma, nu, tau) {
                (-sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)*log(y/nu)*log((nu/mu)^sigma + 1) - 2*sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)*log(nu/mu) + sigma*((y/nu)^sigma + 1)*(-sigma*(nu/mu)^(2*sigma)*((y/nu)^sigma + 1)*log(nu/mu)^2 - (nu/mu)^sigma*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)*log(nu/mu) + ((nu/mu)^sigma + 1)*(-sigma*(y/nu)^sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau)*log(y/nu)*log(nu/mu) + sigma*(y/nu)^sigma*(nu/mu)^sigma*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*(y/nu)^sigma*(nu/mu)^sigma*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log(nu/mu) - 2*sigma*(y/nu)^sigma*(nu/mu)^sigma*log(y/nu)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)^2*log((nu/mu)^sigma + 1) + sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*log(tau)*log(y/nu)*log(nu/y)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu)^2 + sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(tau)*log(nu/mu)*log(nu/y)*log((nu/mu)^sigma + 1) + sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(tau)*log(nu/mu)*log(nu/y) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)^2*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*((y/nu)^sigma*log(y/nu) + ((y/nu)^sigma + 1)*log(nu/y))*log(tau)*log(nu/mu) + (y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) + (y/nu)^sigma*((nu/mu)^sigma + 1)*log(y/nu)*log((nu/mu)^sigma + 1)^2 - (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1)^2 + (nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1)))*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*(-sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) + sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) + sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) - sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) - ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)*log((nu/mu)^sigma + 1))/(sigma^2*((y/nu)^sigma + 1)^2*((nu/mu)^sigma + 1)^2*log((nu/mu)^sigma + 1)^3)
            },
              d2ldddv = function(y, mu, sigma, nu, tau) {
                ((y/nu)^sigma*((nu/mu)^sigma + 1)*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)*log((nu/mu)^sigma + 1) - (nu/mu)^sigma*((y/nu)^sigma + 1)*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2)*log((nu/mu)^sigma + 1) - 2*(nu/mu)^sigma*((y/nu)^sigma + 1)*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) + sigma*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2) + ((y/nu)^sigma + 1)*(-sigma*(nu/mu)^(2*sigma)*((y/nu)^sigma + 1)*log(nu/mu) - sigma*(nu/mu)^sigma*((nu/mu)^sigma + 1)*log(tau)*log(nu/mu) + ((nu/mu)^sigma + 1)*(sigma*(y/nu)^sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau)*log(nu/mu) + sigma*(y/nu)^sigma*(nu/mu)^sigma*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) + sigma*(y/nu)^sigma*(nu/mu)^sigma*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu) - sigma*(y/nu)^sigma*(nu/mu)^sigma*log(y/nu)*log((nu/mu)^sigma + 1) + sigma*(y/nu)^sigma*(nu/mu)^sigma*log(nu/mu)*log((nu/mu)^sigma + 1) - sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log(y/nu)*log((nu/mu)^sigma + 1) - sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) - sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau)*log(nu/mu) + sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(tau)*log(nu/y)*log((nu/mu)^sigma + 1) + sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(tau)*log(nu/y) - sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu)*log((nu/mu)^sigma + 1) - (y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) - (y/nu)^sigma*((nu/mu)^sigma + 1)*log((nu/mu)^sigma + 1)^2 - (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1)^2 + (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1)))*log((nu/mu)^sigma + 1))/(nu*((y/nu)^sigma + 1)^2*((nu/mu)^sigma + 1)^2*log((nu/mu)^sigma + 1)^3)
            },
              d2ldmdd = function(y, mu, sigma, nu, tau) {
                (nu/mu)^sigma*(-sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))*log(nu/mu)*log((nu/mu)^sigma + 1) - 2*sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))*log(nu/mu) + sigma*((nu/mu)^sigma*((y/nu)^sigma + 1)*log(nu/mu) + ((nu/mu)^sigma + 1)*((y/nu)^sigma*log(y/nu) + ((y/nu)^sigma + 1)*log(nu/y))*log(tau))*log((nu/mu)^sigma + 1) + (sigma*log(nu/mu) + 1)*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1))/(mu*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)^2*log((nu/mu)^sigma + 1)^3)
            },
               d2ldv2 = function(y, mu, sigma, nu, tau) {
                sigma*(sigma*(y/nu)^sigma*((nu/mu)^sigma + 1)*(-(y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) - (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) - (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) + sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*((y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) + (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1) - ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) + 2*sigma*(nu/mu)^sigma*((y/nu)^sigma + 1)*((y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) + (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1) - ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1)) + sigma*((y/nu)^sigma + 1)*(-(nu/mu)^(2*sigma)*((y/nu)^sigma + 1) - (nu/mu)^sigma*((nu/mu)^sigma + 1)*log(tau) + ((nu/mu)^sigma + 1)*((y/nu)^sigma*(nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) - (y/nu)^sigma*(nu/mu)^sigma*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) - (y/nu)^sigma*(nu/mu)^sigma*(log(tau) - log((nu/mu)^sigma + 1)) + 2*(y/nu)^sigma*(nu/mu)^sigma*log((nu/mu)^sigma + 1) + (y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) - (y/nu)^sigma*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1) - (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log(tau) - (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1)))*log((nu/mu)^sigma + 1) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*((y/nu)^sigma*((nu/mu)^sigma + 1)*(log(tau) - log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) + (nu/mu)^sigma*(sigma*log(nu/y) + log((y/nu)^sigma + 1))*((y/nu)^sigma + 1)*log(tau) + (nu/mu)^sigma*((y/nu)^sigma + 1)*log((nu/mu)^sigma + 1) - ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*log(tau)*log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1))/(nu^2*((y/nu)^sigma + 1)^2*((nu/mu)^sigma + 1)^2*log((nu/mu)^sigma + 1)^3)
            },
              d2ldmdv = function(y, mu, sigma, nu, tau) {
                sigma^2*(nu/mu)^sigma*(-(nu/mu)^sigma*((y/nu)^sigma + 1)*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) - 2*(nu/mu)^sigma*((y/nu)^sigma + 1)*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1)) + ((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1) + ((nu/mu)^sigma*((y/nu)^sigma + 1) + ((nu/mu)^sigma + 1)*log(tau))*log((nu/mu)^sigma + 1))/(nu*mu*((y/nu)^sigma + 1)*((nu/mu)^sigma + 1)^2*log((nu/mu)^sigma + 1)^3)
            },
            d2ldm2 = function(y, mu, sigma, nu, tau) {
                sigma*(nu/mu)^sigma*(2*sigma*(nu/mu)^sigma*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1)) + sigma*(nu/mu)^sigma*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1) - 1)*log((nu/mu)^sigma + 1) - ((nu/mu)^sigma + 1)*(sigma*((sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1)) + (sigma*log(nu/y) + log((y/nu)^sigma + 1))*log(tau) + log((nu/mu)^sigma + 1))*log((nu/mu)^sigma + 1))/(mu^2*((nu/mu)^sigma + 1)^2*log((nu/mu)^sigma + 1)^3)
            },
          G.dev.incr  = function(y,mu,sigma,nu,tau,...) -2*dRDQ2(y, mu, sigma, nu, tau, log = TRUE), 
                rqres = expression(rqres(pfun="pRDQ2", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)),
           mu.initial = expression(mu <- rep(median(y),length(y))), 
        sigma.initial = expression(sigma <- rep(1, length(y))), 
           nu.initial = expression(nu <- rep(1, length(y))),
          tau.initial = expression(tau <- rep(.5, length(y))),
             mu.valid = function(mu) all(mu > 0), 
          sigma.valid = function(sigma)  all(sigma > 0),
             nu.valid = function(nu) all(nu > 0),
            tau.valid = function(tau) all(tau > 0 & tau < 1),
              y.valid = function(y)  all(y > 0),
                 mean = function(mu, sigma, nu, tau) {
                    a = sigma
                    b = nu
                    p = -log(tau) / log((nu / mu)^sigma + 1)
                    return(ifelse(a > 1, p * b * beta(p+1/a, 1-1/a), Inf))
                  },
             variance = function(mu, sigma, nu, tau) {
                a = sigma
                b = nu
                p = -log(tau) / log((nu / mu)^sigma + 1)
                return(ifelse(a > 2, b^2 * (beta(1-2/a, p+2/a) - beta(1-1/a, p+1/a)^2), Inf))
              }
          ),
                class = c("gamlss.family","family"))
}

dRDQ2 <- function(x, mu = 1, sigma = 1, nu = 1.5, tau = .5, log = FALSE) {
  mu2 <- -log(tau) / log((nu / mu)^sigma + 1)
  ddagum(x = x, scale = nu, shape1.a = sigma, shape2.p = mu2, log = log)
}

pRDQ2 <- function(q, mu = 1, sigma = 1, nu = 1.5, tau = .5, lower.tail = TRUE, log.p = FALSE) {
  mu2 <- -log(tau) / log((nu / mu)^sigma + 1)
  pdagum(q = q, scale = nu, shape1.a = sigma, shape2.p = mu2, lower.tail = lower.tail, log.p = log.p)
}

qRDQ2 <- function(p, mu = 1, sigma = 1, nu = 1.5, tau = .5, lower.tail = TRUE, log.p = FALSE) {
  mu2 <- -log(tau) / log((nu / mu)^sigma + 1)
  qdagum(p = p, scale = nu, shape1.a = sigma, shape2.p = mu2, lower.tail = lower.tail, log.p = log.p)
}

rRDQ2 <- function(n, mu = 1, sigma = 1, nu = 1.5, tau = .5) {
  mu2 <- -log(tau) / log((nu / mu)^sigma + 1)
  rdagum(n = n, scale = nu, shape1.a = sigma, shape2.p = mu2)
}
