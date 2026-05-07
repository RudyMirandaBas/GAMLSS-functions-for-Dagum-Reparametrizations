library(gamlss)
library(VGAM)

logshiftto1 <- function() {
linkfun = function(mu) log(mu - 1)
linkinv = function(eta) 1 + pmax(.Machine$double.eps, exp(eta))
mu.eta = function(eta) pmax(.Machine$double.eps, exp(eta))
valideta = function(eta) TRUE   
#make.link.gamlss("logshiftto1")
link <- "logshiftto1"
structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
        valideta = valideta, name = link), class = "link-gamlss")
}


RD <- function(mu.link="log", sigma.link=logshiftto1(), nu.link="log") {
    mstats <- checklink("mu.link", "Reparametrized.Dagum1", substitute(mu.link), c("inverse", "log", "identity"))# dummy
    dstats <- checklink("sigma.link", "Reparametrized.Dagum1", substitute(sigma.link), c("logshiftto1", "inverse", "log", "identity"))
    vstats <- checklink("nu.link", "Reparametrized.Dagum1", substitute(nu.link), c("inverse", "log", "identity"))

    structure(
          list(family = c("RD", "Reparametrized.Dagum1"),
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
                 dldm = function(y, mu, sigma, nu) {
                  sigma*(-nu + (sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu +
	 1)))^sigma)/(mu*((sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + 1))
                },
               d2ldm2 = function(y, mu, sigma, nu) {
                -sigma*(sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma - nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma - nu + (sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma) + (sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu +
	 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma)/(mu^2*((sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + 1)^2)
              },
                 dldv = function(y, mu, sigma, nu) {
                  -mu*sigma*(nu*y*gamma(1 - 1/sigma)*gamma(nu + 1/sigma)/(mu*gamma(nu + 1)))^sigma*(nu + 1)*(-nu*y*gamma(1 - 1/sigma)*gamma(nu + 1/sigma)*digamma(nu + 1)/(mu*gamma(nu + 1)) + nu*y*gamma(1 - 1/sigma)*gamma(nu + 1/sigma)*digamma(nu + 1/sigma)/(mu*gamma(nu + 1)) + y*gamma(1 - 1/sigma)*gamma(nu + 1/sigma)/(mu*gamma(nu + 1)))*gamma(nu + 1)/(nu*y*((nu*y*gamma(1 - 1/sigma)*gamma(nu + 1/sigma)/(mu*gamma(nu + 1)))^sigma + 1)*gamma(1 - 1/sigma)*gamma(nu + 1/sigma)) - sigma*nu*(gamma(1 - 1/sigma)*gamma(nu + 1/sigma)*digamma(nu + 1)/gamma(nu + 1) - gamma(1 - 1/sigma)*gamma(nu + 1/sigma)*digamma(nu + 1/sigma)/gamma(nu + 1))*gamma(nu + 1)/(gamma(1 - 1/sigma)*gamma(nu + 1/sigma)) - sigma*log(mu) + sigma*log(nu) + sigma*log(y) + sigma*log(gamma(1 - 1/sigma)*gamma(nu + 1/sigma)/gamma(nu + 1)) - log((nu*y*gamma(1 - 1/sigma)*gamma(nu + 1/sigma)/(mu*gamma(nu + 1)))^sigma + 1) + (sigma*nu + 1)/nu
                },
               d2ldv2 = function(y, mu, sigma, nu) {
                (-sigma^2*nu^3*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma)^2 + 2*sigma^2*nu^3*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma)*digamma(nu + 1) - sigma^2*nu^3*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma(nu + 1)^2 - sigma^2*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma)^2 + 2*sigma^2*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma)*digamma(nu + 1) - 2*sigma^2*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma) - sigma^2*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma(nu + 1)^2 + 2*sigma^2*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma(nu + 1) - 2*sigma^2*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma) + 2*sigma^2*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma(nu + 1) - sigma^2*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma - sigma^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + sigma*nu^3*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*trigamma((sigma*nu + 1)/sigma) - sigma*nu^3*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*trigamma(nu + 1) + sigma*nu^3*trigamma((sigma*nu + 1)/sigma) - sigma*nu^3*trigamma(nu + 1) - sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma)*trigamma((sigma*nu + 1)/sigma) + sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma)*trigamma(nu + 1) + 2*sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma) - 2*sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma(nu + 1) - sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*trigamma((sigma*nu + 1)/sigma) + sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*trigamma(nu + 1) + 2*sigma*nu^2*digamma((sigma*nu + 1)/sigma) - 2*sigma*nu^2*digamma(nu + 1) + sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + 
	sigma*nu + sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma) + sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma - (sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma) - 2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma - 1)/(nu^2*((sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + 1)^2)
              },
                 dldd = function(y, mu, sigma, nu) {
                  sigma*nu*(gamma(1 - 1/sigma)*gamma(nu + 1/sigma)*digamma(1 - 1/sigma)/(sigma^2*gamma(nu + 1)) - gamma(1 - 1/sigma)*gamma(nu + 1/sigma)*digamma(nu + 1/sigma)/(sigma^2*gamma(nu + 1)))*gamma(nu + 1)/(gamma(1 - 1/sigma)*gamma(nu + 1/sigma)) - nu*log(mu) + nu*log(nu) + nu*log(y) + nu*log(gamma(1 - 1/sigma)*gamma(nu + 1/sigma)/gamma(nu + 1)) - (nu*y*gamma(1 - 1/sigma)*gamma(nu + 1/sigma)/(mu*gamma(nu + 1)))^sigma*(nu + 1)*(mu*sigma*(nu*y*gamma(1 - 1/sigma)*gamma(nu + 1/sigma)*digamma(1 - 1/sigma)/(mu*sigma^2*gamma(nu + 1)) - nu*y*gamma(1 - 1/sigma)*gamma(nu + 1/sigma)*digamma(nu + 1/sigma)/(mu*sigma^2*gamma(nu + 1)))*gamma(nu + 1)/(nu*y*gamma(1 - 1/sigma)*gamma(nu + 1/sigma)) + log(nu*y*gamma(1 - 1/sigma)*gamma(nu + 1/sigma)/(mu*gamma(nu + 1))))/((nu*y*gamma(1 - 1/sigma)*gamma(nu + 1/sigma)/(mu*gamma(nu + 1)))^sigma + 1) + 1/sigma
                },
               d2ldd2 = function(y, mu, sigma, nu) {
                (-sigma^3*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu)))^2 - sigma^3*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu)))^2 - 2*sigma^2*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu)))*digamma((sigma - 1)/sigma) + 2*sigma^2*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu)))*digamma((sigma*nu + 1)/sigma) - 2*sigma^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu)))*digamma((sigma - 1)/sigma) + 2*sigma^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu)))*digamma((sigma*nu + 1)/sigma) - sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma)^2 + 2*sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma)*digamma((sigma*nu + 1)/sigma) - sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma)^2 - sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma) - sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma)^2 + 2*sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma)*digamma((sigma*nu + 1)/sigma) - sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma)^2 - 2*sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma - sigma + nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*trigamma((sigma - 1)/sigma) + nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*trigamma((sigma*nu + 1)/sigma) + nu*trigamma((sigma - 1)/sigma) + nu*trigamma((sigma*nu + 1)/sigma) - (sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma)*trigamma((sigma - 1)/sigma) - (sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma)*trigamma((sigma*nu + 1)/sigma) - (sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*trigamma((sigma - 
	1)/sigma) - (sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*trigamma((sigma*nu + 1)/sigma))/(sigma^3*((sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + 1)^2)
              },
              d2ldmdv = function(y, mu, sigma, nu) {
                -sigma*(-sigma^2*nu^3*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + sigma + 1)/sigma) + sigma^2*nu^3*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma(nu + 1) - sigma^2*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + sigma + 1)/sigma) + sigma^2*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma(nu + 1) - sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + sigma + 1)/sigma) + sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma(nu + 1) + sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + sigma*nu^2 - sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + sigma + 1)/sigma) + sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma(nu + 1) - sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma - sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu +
	 1)))^sigma + nu)/(mu*nu*(sigma*nu + 1)*((sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + 1)^2)
              },
              d2ldmdd = function(y, mu, sigma, nu) {
                (sigma^2*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu))) + sigma^2*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu))) + sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma) - sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + sigma + 1)/sigma) - sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma - sigma*nu^2 + sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma) + sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu))) + sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma) - sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + sigma + 1)/sigma) + 2*sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu))) + sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma) - nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + sigma + 1)/sigma) - nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma - nu + (sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma) + (sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma) - (sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + sigma + 1)/sigma) + (sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 
	1)))^sigma)/(mu*(sigma*nu + 1)*((sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + 1)^2)
              },
              d2ldddv = function(y, mu, sigma, nu) {
                (-sigma^2*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu)))*digamma((sigma*nu + 1)/sigma) + sigma^2*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu)))*digamma(nu + 1) - sigma^2*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu)))*digamma((sigma*nu + 1)/sigma) + sigma^2*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu)))*digamma(nu + 1) - sigma^2*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu))) - sigma^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu))) - sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma)*digamma((sigma*nu + 1)/sigma) + sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma)*digamma(nu + 1) + sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma)^2 - sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma)*digamma(nu + 1) + sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma) - sigma*nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma(nu + 1) + sigma*nu^2*digamma((sigma*nu + 1)/sigma) - sigma*nu^2*digamma(nu + 1) - sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma)*log(mu) + sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma)*log(nu) + sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma)*log(y) + sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma)*log(sigma*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/((sigma*nu + 1)*gamma(nu + 1))) - sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma)*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu))) - sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma)*digamma((sigma*nu +
	 1)/sigma) + sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma)*digamma(nu + 1) - 2*sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(mu) + 2*sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(nu) + 2*sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(y) + 2*sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/((sigma*nu + 1)*gamma(nu + 1))) - sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*log(sigma*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu))) - sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma)*digamma((sigma*nu + 1)/sigma) + sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma)*digamma(nu + 1) - sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma) + sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma)^2 - sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma)*digamma(nu + 1) + sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma(nu + 1) + sigma*nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma - sigma*nu*log(mu) + sigma*nu*log(nu) + sigma*nu*log(y) + sigma*nu*log(sigma*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/((sigma*nu + 1)*gamma(nu + 1))) + sigma*nu - sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma) - sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma) + sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma) - sigma*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma - nu^2*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*trigamma((sigma*nu + 1)/sigma) - nu^2*trigamma((sigma*nu + 1)/sigma) + nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^(2*sigma)*trigamma((sigma*nu + 1)/sigma) + nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma - 1)/sigma) - nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*digamma((sigma*nu + 1)/sigma) + nu*(sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 + 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma*trigamma((sigma*nu + 1)/sigma) + nu*digamma((sigma - 1)/sigma) - nu*digamma((sigma*nu + 1)/sigma))/(sigma*nu*((sigma*nu*y*gamma(1 - 1/sigma)*gamma(nu + 1 +
	 1/sigma)/(mu*(sigma*nu + 1)*gamma(nu + 1)))^sigma + 1)^2)
              },
          G.dev.incr  = function(y,mu,sigma,nu,...) -2*dRD(y, mu, sigma, nu, log = TRUE), 
                rqres = expression(rqres(pfun="pRD", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
           mu.initial = expression(mu <- rep(1,length(y))),
        sigma.initial = expression(sigma <- rep(1.5, length(y))),
           nu.initial = expression(nu <- rep(1, length(y))),
             mu.valid = function(mu) all(mu > 0),
          sigma.valid = function(sigma) all(sigma > 1),
             nu.valid = function(nu)  all(nu > 0),
              y.valid = function(y)  all(y > 0),
                 mean = function(mu) ifelse(sigma > 1, mu, Inf),
             variance = function(mu, sigma, nu) {
                a = sigma
                b = exp(log(mu) - log(nu) - lbeta(nu+1/sigma, 1-1/sigma))
                p = nu
                return(ifelse(a > 2, b^2 * (beta(1-2/a, p+2/a) - beta(1-1/a, p+1/a)^2), Inf))
              }
          ),
                class = c("gamlss.family","family"))
}

dRD <- function(x, mu = 1, sigma = 1.5, nu = 1, log = FALSE) {
  mu2 <- exp(log(mu) - log(nu) - lbeta(nu+1/sigma, 1-1/sigma))
  fx <- ddagum(x = x, scale = mu2, shape1.a = sigma, shape2.p = nu, log = log)
  ifelse(x > 0, fx, 0)
}

pRD <- function(q, mu = 1, sigma = 1.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  mu2 <- exp(log(mu) - log(nu) - lbeta(nu+1/sigma, 1-1/sigma))
  pdagum(q = q, scale = mu2, shape1.a = sigma, shape2.p = nu, lower.tail = lower.tail, log.p = log.p)
}

qRD <- function(p, mu = 1, sigma = 1.5, nu = 1, lower.tail = TRUE, log.p = FALSE) {
  mu2 <- exp(log(mu) - log(nu) - lbeta(nu+1/sigma, 1-1/sigma))
  qdagum(p = p, scale = mu2, shape1.a = sigma, shape2.p = nu, lower.tail = lower.tail, log.p = log.p)
}

rRD <- function(n, mu = 1, sigma = 1.5, nu = 1) {
  mu2 <- exp(log(mu) - log(nu) - lbeta(nu+1/sigma, 1-1/sigma))
  rdagum(n = n, scale = mu2, shape1.a = sigma, shape2.p = nu)
}

#dat <- rRD(1000, 5, 4, 1.8)
#gamlss(dat~1, family=RD)
#gamlss(dat~1, family=WEI)
#
#l <- function(theta, x) {
#  mu=theta[1]; nu=theta[2]; sigma=theta[3]
#  # Define B
#  B = beta(nu+1/sigma, 1-1/sigma)
#
#  # Define L
#  (1+nu*sigma)*log(nu) + log(sigma) + (nu*sigma-1)*log(x) - nu*sigma*log(mu) + nu*sigma*log(B) - (nu+1)*log(1+(x*nu*B/mu)^sigma)
#}
#
#optim(par = c(5,1,1.5), fn = function(x) -sum(dRD(x = dat, mu = x[1], nu = x[2], sigma = x[3], log = TRUE)), control=list(maxit=10000))
