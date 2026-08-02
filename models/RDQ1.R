library(gamlss)
library(VGAM)

RDQ1 <- function (mu.link="log", sigma.link="log", nu.link="log", tau.link="identity") {
    mstats <- checklink("mu.link", "Reparametrized.Dagum1", substitute(mu.link), c("inverse", "log", "identity"))# dummy
    dstats <- checklink("sigma.link", "Reparametrized.Dagum1", substitute(sigma.link), c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "Reparametrized.Dagum1", substitute(nu.link), c("inverse", "log", "identity"))
    tstats <- checklink("tau.link", "Reparametrized.Dagum1", substitute(tau.link), c("inverse", "log", "identity"))


    structure(
          list(family = c("RDQ1", "Reparametrized.Dagum1"),
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
               dldm = function(y, mu, sigma, nu, tau) {
                    return(sigma*(nu*tau^(1/nu) - nu + tau^(1/nu)*(y/mu)^sigma)/(mu*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)))
                },
                 dldd = function(y, mu, sigma, nu, tau) {
                    return((nu*sigma*(-log(mu) + log(y))*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1) - sigma*tau^(1/nu)*(y/mu)^sigma*(nu + 1)*log(y/mu) + tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)/(sigma*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)))
                },
                 dldv = function(y, mu, sigma, nu, tau) {
                    return((nu^2*(tau^(1/nu) - 1)*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)*(-sigma*log(mu) + sigma*log(y) - log((1 - tau^(1/nu))/tau^(1/nu)) - log((-tau^(1/nu)*(y/mu)^sigma + tau^(1/nu) - 1)/(tau^(1/nu) - 1))) + nu*(tau^(1/nu) - 1)*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1) + nu*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)*log(tau) - tau^(1/nu)*(y/mu)^sigma*(nu + 1)*log(tau))/(nu^2*(tau^(1/nu) - 1)*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)))
                },
               d2ldm2 = function(y, mu, sigma, nu, tau) {
                return(sigma*(sigma*tau^(1/nu)*(y/mu)^sigma*(nu*tau^(1/nu) - nu + tau^(1/nu)*(y/mu)^sigma) + (tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)*(-nu*tau^(1/nu) + nu - sigma*tau^(1/nu)*(y/mu)^sigma - tau^(1/nu)*(y/mu)^sigma))/(mu^2*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)^2))
            },
              d2ldmdd = function(y, mu, sigma, nu, tau) {
                return((-sigma*tau^(1/nu)*(y/mu)^sigma*(nu*tau^(1/nu) - nu + tau^(1/nu)*(y/mu)^sigma)*log(y/mu) + (tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)*(nu*tau^(1/nu) - nu + sigma*tau^(1/nu)*(y/mu)^sigma*log(y/mu) + tau^(1/nu)*(y/mu)^sigma))/(mu*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)^2))
            },
              d2ldmdv = function(y, mu, sigma, nu, tau) {
                return(sigma*(tau^(1/nu)*((y/mu)^sigma - 1)*(nu*tau^(1/nu) - nu + tau^(1/nu)*(y/mu)^sigma)*log(tau) - (nu^2*(1 - tau^(1/nu)) + nu*tau^(1/nu)*log(tau) + tau^(1/nu)*(y/mu)^sigma*log(tau))*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1))/(mu*nu^2*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)^2))
            },
               d2ldd2 = function(y, mu, sigma, nu, tau) {
                return((nu*sigma^2*tau^(2/nu)*(y/mu)^sigma*log(y/mu)^2 - nu*sigma^2*tau^(1/nu)*(y/mu)^sigma*log(y/mu)^2 + sigma^2*tau^(2/nu)*(y/mu)^sigma*log(y/mu)^2 - sigma^2*tau^(1/nu)*(y/mu)^sigma*log(y/mu)^2 - tau^(2/nu)*(y/mu)^(2*sigma) + 2*tau^(2/nu)*(y/mu)^sigma - tau^(2/nu) - 2*tau^(1/nu)*(y/mu)^sigma + 2*tau^(1/nu) - 1)/(sigma^2*(tau^(2/nu)*(y/mu)^(2*sigma) - 2*tau^(2/nu)*(y/mu)^sigma + tau^(2/nu) + 2*tau^(1/nu)*(y/mu)^sigma - 2*tau^(1/nu) + 1)))
            },
              d2ldddv = function(y, mu, sigma, nu, tau) {
                return((-tau^(1/nu)*((y/mu)^sigma - 1)*(nu*sigma*(log(mu) - log(y))*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1) + sigma*tau^(1/nu)*(y/mu)^sigma*(nu + 1)*log(y/mu) - tau^(1/nu)*(y/mu)^sigma + tau^(1/nu) - 1)*log(tau) + (tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)*(-nu^2*sigma*(tau^(1/nu)*(y/mu)^sigma*log(y/mu) + (log(mu) - log(y))*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)) + nu*sigma*tau^(1/nu)*((y/mu)^sigma - 1)*(log(mu) - log(y))*log(tau) + tau^(1/nu)*(sigma*(y/mu)^sigma*(nu + 1)*log(y/mu) - (y/mu)^sigma + 1)*log(tau)))/(nu^2*sigma*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)^2))
            },
            d2ldv2 = function(y, mu, sigma, nu, tau) {
                return((-nu^2*tau^(4/nu)*(y/mu)^(2*sigma) + 2*nu^2*tau^(4/nu)*(y/mu)^sigma - nu^2*tau^(4/nu) + 2*nu^2*tau^(3/nu)*(y/mu)^(2*sigma) - 6*nu^2*tau^(3/nu)*(y/mu)^sigma + 4*nu^2*tau^(3/nu) - nu^2*tau^(2/nu)*(y/mu)^(2*sigma) + 6*nu^2*tau^(2/nu)*(y/mu)^sigma - 6*nu^2*tau^(2/nu) - 2*nu^2*tau^(1/nu)*(y/mu)^sigma + 4*nu^2*tau^(1/nu) - nu^2 + 2*nu*tau^(3/nu)*(y/mu)^(2*sigma)*log(tau) - nu*tau^(3/nu)*(y/mu)^sigma*log(tau)^2 - 2*nu*tau^(3/nu)*(y/mu)^sigma*log(tau) + nu*tau^(3/nu)*log(tau)^2 - 2*nu*tau^(2/nu)*(y/mu)^(2*sigma)*log(tau) + 2*nu*tau^(2/nu)*(y/mu)^sigma*log(tau)^2 + 4*nu*tau^(2/nu)*(y/mu)^sigma*log(tau) - 2*nu*tau^(2/nu)*log(tau)^2 - nu*tau^(1/nu)*(y/mu)^sigma*log(tau)^2 - 2*nu*tau^(1/nu)*(y/mu)^sigma*log(tau) + nu*tau^(1/nu)*log(tau)^2 - tau^(3/nu)*(y/mu)^(2*sigma)*log(tau)^2 + tau^(3/nu)*(y/mu)^sigma*log(tau)^2 - tau^(1/nu)*(y/mu)^sigma*log(tau)^2)/(nu^4*(tau^(4/nu)*(y/mu)^(2*sigma) - 2*tau^(4/nu)*(y/mu)^sigma + tau^(4/nu) - 2*tau^(3/nu)*(y/mu)^(2*sigma) + 6*tau^(3/nu)*(y/mu)^sigma - 4*tau^(3/nu) + tau^(2/nu)*(y/mu)^(2*sigma) - 6*tau^(2/nu)*(y/mu)^sigma + 6*tau^(2/nu) + 2*tau^(1/nu)*(y/mu)^sigma - 4*tau^(1/nu) + 1)))
            },
          G.dev.incr  = function(y,mu,sigma,nu,tau,...) -2*dRDQ1(y, mu, sigma, nu, tau, log = TRUE), 
                rqres = expression(rqres(pfun="pRDQ1", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)),
           mu.initial = expression(mu <- rep(median(y),length(y))), 
        sigma.initial = expression(sigma <- rep(1, length(y))), 
           nu.initial = expression(nu <- rep(1, length(y))),
          tau.initial = expression(tau <- rep(.5, length(y))),
             mu.valid = function(mu) all(mu > 0), 
          sigma.valid = function(sigma)  all(sigma > 0),
             nu.valid = function(nu) all(nu > 0),
            tau.valid = function(tau) all(tau > 0),
              y.valid = function(y)  all(y > 0),
                 mean = function(mu, sigma, nu, tau) {
                    a = sigma
                    b = mu * (tau^(-1/nu)-1)^(1/sigma)
                    p = nu
                    return(ifelse(a > 1, p * b * beta(p+1/a, 1-1/a), Inf))
                  },
             variance = function(mu, sigma, nu, tau) {
                a = sigma
                b = mu * (tau^(-1/nu)-1)^(1/sigma)
                p = nu
                return(ifelse(a > 2, b^2 * (beta(1-2/a, p+2/a) - beta(1-1/a, p+1/a)^2), Inf))
              }
          ),
                class = c("gamlss.family","family"))
}

dRDQ1 <- function(x, mu = 1, sigma = 1, nu = 1.5, tau = .5, log = FALSE) {
  mu2 <- mu * (tau^(-1 / nu) - 1)^(1 / sigma)
  ddagum(x = x, scale = mu2, shape1.a = sigma, shape2.p = nu, log = log)
}

pRDQ1 <- function(q, mu = 1, sigma = 1, nu = 1.5, tau = .5, lower.tail = TRUE, log.p = FALSE) {
  mu2 <- mu * (tau^(-1 / nu) - 1)^(1 / sigma)
  pdagum(q = q, scale = mu2, shape1.a = sigma, shape2.p = nu, lower.tail = lower.tail, log.p = log.p)
}

qRDQ1 <- function(p, mu = 1, sigma = 1, nu = 1.5, tau = .5, lower.tail = TRUE, log.p = FALSE) {
  mu2 <- mu * (tau^(-1 / nu) - 1)^(1 / sigma)
  qdagum(p = p, scale = mu2, shape1.a = sigma, shape2.p = nu, lower.tail = lower.tail, log.p = log.p)
}

rRDQ1 <- function(n, mu = 1, sigma = 1, nu = 1.5, tau = .5) {
  mu2 <- mu * (tau^(-1 / nu) - 1)^(1 / sigma)
  rdagum(n = n, scale = mu2, shape1.a = sigma, shape2.p = nu)
}
