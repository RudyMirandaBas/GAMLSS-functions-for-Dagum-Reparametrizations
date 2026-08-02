library(gamlss)
library(VGAM)

set.seed(123)

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


RD <- function(mu.link="log", sigma.link="log", nu.link=logshiftto1()) {
    mstats <- checklink("mu.link", "Reparametrized.Dagum1", substitute(mu.link), c("inverse", "log", "identity"))# dummy
    dstats <- checklink("sigma.link", "Reparametrized.Dagum1", substitute(sigma.link), c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "Reparametrized.Dagum1", substitute(nu.link), c("logshiftto1", "inverse", "log", "identity"))

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
                  nu*(-sigma + (nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma +
	 1)))^nu)/(mu*((nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + 1))
                },
               d2ldm2 = function(y, mu, sigma, nu) {
                -nu*(nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu - sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu - sigma + (nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu) + (nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma +
	 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu)/(mu^2*((nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + 1)^2)
              },
                 dldd = function(y, mu, sigma, nu) {
                  -mu*nu*(sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1/nu)/(mu*gamma(sigma + 1)))^nu*(sigma + 1)*(-sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1/nu)*digamma(sigma + 1)/(mu*gamma(sigma + 1)) + sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1/nu)*digamma(sigma + 1/nu)/(mu*gamma(sigma + 1)) + y*gamma(1 - 1/nu)*gamma(sigma + 1/nu)/(mu*gamma(sigma + 1)))*gamma(sigma + 1)/(sigma*y*((sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1/nu)/(mu*gamma(sigma + 1)))^nu + 1)*gamma(1 - 1/nu)*gamma(sigma + 1/nu)) - nu*sigma*(gamma(1 - 1/nu)*gamma(sigma + 1/nu)*digamma(sigma + 1)/gamma(sigma + 1) - gamma(1 - 1/nu)*gamma(sigma + 1/nu)*digamma(sigma + 1/nu)/gamma(sigma + 1))*gamma(sigma + 1)/(gamma(1 - 1/nu)*gamma(sigma + 1/nu)) - nu*log(mu) + nu*log(sigma) + nu*log(y) + nu*log(gamma(1 - 1/nu)*gamma(sigma + 1/nu)/gamma(sigma + 1)) - log((sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1/nu)/(mu*gamma(sigma + 1)))^nu + 1) + (nu*sigma + 1)/sigma
                },
               d2ldd2 = function(y, mu, sigma, nu) {
                (-nu^2*sigma^3*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu)^2 + 2*nu^2*sigma^3*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu)*digamma(sigma + 1) - nu^2*sigma^3*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma(sigma + 1)^2 - nu^2*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu)^2 + 2*nu^2*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu)*digamma(sigma + 1) - 2*nu^2*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu) - nu^2*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma(sigma + 1)^2 + 2*nu^2*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma(sigma + 1) - 2*nu^2*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu) + 2*nu^2*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma(sigma + 1) - nu^2*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu - nu^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + nu*sigma^3*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*trigamma((nu*sigma + 1)/nu) - nu*sigma^3*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*trigamma(sigma + 1) + nu*sigma^3*trigamma((nu*sigma + 1)/nu) - nu*sigma^3*trigamma(sigma + 1) - nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu)*trigamma((nu*sigma + 1)/nu) + nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu)*trigamma(sigma + 1) + 2*nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu) - 2*nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma(sigma + 1) - nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*trigamma((nu*sigma + 1)/nu) + nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*trigamma(sigma + 1) + 2*nu*sigma^2*digamma((nu*sigma + 1)/nu) - 2*nu*sigma^2*digamma(sigma + 1) + nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + 
	nu*sigma + nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu) + nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu - (nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu) - 2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu - 1)/(sigma^2*((nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + 1)^2)
              },
                 dldv = function(y, mu, sigma, nu) {
                  nu*sigma*(gamma(1 - 1/nu)*gamma(sigma + 1/nu)*digamma(1 - 1/nu)/(nu^2*gamma(sigma + 1)) - gamma(1 - 1/nu)*gamma(sigma + 1/nu)*digamma(sigma + 1/nu)/(nu^2*gamma(sigma + 1)))*gamma(sigma + 1)/(gamma(1 - 1/nu)*gamma(sigma + 1/nu)) - sigma*log(mu) + sigma*log(sigma) + sigma*log(y) + sigma*log(gamma(1 - 1/nu)*gamma(sigma + 1/nu)/gamma(sigma + 1)) - (sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1/nu)/(mu*gamma(sigma + 1)))^nu*(sigma + 1)*(mu*nu*(sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1/nu)*digamma(1 - 1/nu)/(mu*nu^2*gamma(sigma + 1)) - sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1/nu)*digamma(sigma + 1/nu)/(mu*nu^2*gamma(sigma + 1)))*gamma(sigma + 1)/(sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1/nu)) + log(sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1/nu)/(mu*gamma(sigma + 1))))/((sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1/nu)/(mu*gamma(sigma + 1)))^nu + 1) + 1/nu
                },
               d2ldv2 = function(y, mu, sigma, nu) {
                (-nu^3*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma)))^2 - nu^3*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma)))^2 - 2*nu^2*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma)))*digamma((nu - 1)/nu) + 2*nu^2*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma)))*digamma((nu*sigma + 1)/nu) - 2*nu^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma)))*digamma((nu - 1)/nu) + 2*nu^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma)))*digamma((nu*sigma + 1)/nu) - nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu)^2 + 2*nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu)*digamma((nu*sigma + 1)/nu) - nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu)^2 - nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu) - nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu)^2 + 2*nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu)*digamma((nu*sigma + 1)/nu) - nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu)^2 - 2*nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu - nu + sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*trigamma((nu - 1)/nu) + sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*trigamma((nu*sigma + 1)/nu) + sigma*trigamma((nu - 1)/nu) + sigma*trigamma((nu*sigma + 1)/nu) - (nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu)*trigamma((nu - 1)/nu) - (nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu)*trigamma((nu*sigma + 1)/nu) - (nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*trigamma((nu - 
	1)/nu) - (nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*trigamma((nu*sigma + 1)/nu))/(nu^3*((nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + 1)^2)
              },
              d2ldmdd = function(y, mu, sigma, nu) {
                -nu*(-nu^2*sigma^3*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + nu + 1)/nu) + nu^2*sigma^3*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma(sigma + 1) - nu^2*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + nu + 1)/nu) + nu^2*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma(sigma + 1) - nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + nu + 1)/nu) + nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma(sigma + 1) + nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + nu*sigma^2 - nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + nu + 1)/nu) + nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma(sigma + 1) - nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu - nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma +
	 1)))^nu + sigma)/(mu*sigma*(nu*sigma + 1)*((nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + 1)^2)
              },
              d2ldmdv = function(y, mu, sigma, nu) {
                (nu^2*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma))) + nu^2*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma))) + nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu) - nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + nu + 1)/nu) - nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu - nu*sigma^2 + nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu) + nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma))) + nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu) - nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + nu + 1)/nu) + 2*nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma))) + nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu) - sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + nu + 1)/nu) - sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu - sigma + (nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu) + (nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu) - (nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + nu + 1)/nu) + (nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 
	1)))^nu)/(mu*(nu*sigma + 1)*((nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + 1)^2)
              },
              d2ldddv = function(y, mu, sigma, nu) {
                (-nu^2*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma)))*digamma((nu*sigma + 1)/nu) + nu^2*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma)))*digamma(sigma + 1) - nu^2*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma)))*digamma((nu*sigma + 1)/nu) + nu^2*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma)))*digamma(sigma + 1) - nu^2*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma))) - nu^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma))) - nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu)*digamma((nu*sigma + 1)/nu) + nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu)*digamma(sigma + 1) + nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu)^2 - nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu)*digamma(sigma + 1) + nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu) - nu*sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma(sigma + 1) + nu*sigma^2*digamma((nu*sigma + 1)/nu) - nu*sigma^2*digamma(sigma + 1) - nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu)*log(mu) + nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu)*log(sigma) + nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu)*log(y) + nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu)*log(nu*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/((nu*sigma + 1)*gamma(sigma + 1))) - nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu)*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma))) - nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu)*digamma((nu*sigma +
	 1)/nu) + nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu)*digamma(sigma + 1) - 2*nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(mu) + 2*nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(sigma) + 2*nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(y) + 2*nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/((nu*sigma + 1)*gamma(sigma + 1))) - nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*log(nu*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma))) - nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu)*digamma((nu*sigma + 1)/nu) + nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu)*digamma(sigma + 1) - nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu) + nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu)^2 - nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu)*digamma(sigma + 1) + nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma(sigma + 1) + nu*sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu - nu*sigma*log(mu) + nu*sigma*log(sigma) + nu*sigma*log(y) + nu*sigma*log(nu*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/((nu*sigma + 1)*gamma(sigma + 1))) + nu*sigma - nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu) - nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu) + nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu) - nu*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu - sigma^2*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*trigamma((nu*sigma + 1)/nu) - sigma^2*trigamma((nu*sigma + 1)/nu) + sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^(2*nu)*trigamma((nu*sigma + 1)/nu) + sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu - 1)/nu) - sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*digamma((nu*sigma + 1)/nu) + sigma*(nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 + 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu*trigamma((nu*sigma + 1)/nu) + sigma*digamma((nu - 1)/nu) - sigma*digamma((nu*sigma + 1)/nu))/(nu*sigma*((nu*sigma*y*gamma(1 - 1/nu)*gamma(sigma + 1 +
	 1/nu)/(mu*(nu*sigma + 1)*gamma(sigma + 1)))^nu + 1)^2)
              },
          G.dev.incr  = function(y,mu,sigma,nu,...) -2*dRD(y, mu, sigma, nu, log = TRUE), 
                rqres = expression(rqres(pfun="pRD", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
           mu.initial = expression(mu <- rep(1,length(y))),
        sigma.initial = expression(sigma <- rep(1, length(y))),
           nu.initial = expression(nu <- rep(1.5, length(y))),
             mu.valid = function(mu) all(mu > 0),
          sigma.valid = function(sigma)  all(sigma > 0),
             nu.valid = function(nu) all(nu > 1),
              y.valid = function(y)  all(y > 0),
                 mean = function(mu) ifelse(nu > 1, mu, Inf),
             variance = function(mu, sigma, nu) {
                a = nu
                b = mu/(sigma*beta(sigma+1/nu, 1-1/nu))
                p = sigma
                return(ifelse(a > 2, b^2 * (beta(1-2/a, p+2/a) - beta(1-1/a, p+1/a)^2), Inf))
              }
          ),
                class = c("gamlss.family","family"))
}

dRD <- function(x, mu = 1, sigma = 1, nu = 1.5, log = FALSE) {
  mu2 <- mu/(sigma*beta(sigma+1/nu, 1-1/nu))
  fx <- ddagum(x = x, scale = mu2, shape1.a = nu, shape2.p = sigma, log = log)
  ifelse(x > 0, fx, 0)
}

pRD <- function(q, mu = 1, sigma = 1, nu = 1.5, lower.tail = TRUE, log.p = FALSE) {
  mu2 <- mu/(sigma*beta(sigma+1/nu, 1-1/nu))
  pdagum(q = q, scale = mu2, shape1.a = nu, shape2.p = sigma, lower.tail = lower.tail, log.p = log.p)
}

qRD <- function(p, mu = 1, sigma = 1, nu = 1.5, lower.tail = TRUE, log.p = FALSE) {
  mu2 <- mu/(sigma*beta(sigma+1/nu, 1-1/nu))
  qdagum(p = p, scale = mu2, shape1.a = nu, shape2.p = sigma, lower.tail = lower.tail, log.p = log.p)
}

rRD <- function(n, mu = 1, sigma = 1, nu = 1.5) {
  mu2 <- mu/(sigma*beta(sigma+1/nu, 1-1/nu))
  rdagum(n = n, scale = mu2, shape1.a = nu, shape2.p = sigma)
}

#######################################
############# Monte Carlo #############
#######################################

library(parallel)

MCllikelihood_estimation <- function(n = 1000, mu_ = c(1, .5, .2), sigma_ = c(.5, .4), nu_ = c(-.2)) {

  z1 <- as.numeric(scale(runif(n)))
  z2 <- as.numeric(scale(rnorm(n)))
  
  X1 <- model.matrix(~z1+z2)                  #matriz de diseño de mu
  X2 <- model.matrix(~z1)                     #matriz de diseño de sigma
  X3 <- matrix(1, nrow = n)                     #matriz de diseño de nu
  colnames(X3) <- "(Intercept)"
  
  mu.true <- as.vector(exp(X1 %*% mu_))     # log-link
  sigma.true <- as.vector(exp(X2 %*% sigma_))  # log-link
  nu.true <- as.vector(1 + exp(X3 %*% nu_)) # log-link-desplazado (nu > 1)
  y <- rRD(n, mu.true, sigma.true, nu.true) ##simulando valores
  
  # GAMLSS
  aux <- gamlss(y ~ X1[, -1, drop=FALSE],
                sigma.fo = ~X2[, -1, drop = FALSE],
                family = RD, method = RS(1000),
                control = gamlss.control(trace = FALSE))
  capture.output(res.gamlss <- summary(aux)[, 1:2], file = nullfile())

  res <- cbind(c(mu_, sigma_, nu_),                              # REAL
               res.gamlss)                                       # GAMLSS

  colnames(res) <- c("true",
                     "est gamlss", "se gamlss")
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

prob <- function(i, theta, seeds, casos.n) {
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
                                            nu_ = theta[6]),
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
system.time(testing.p1 <- mclapply(1:4000, prob,
                                   seeds = seeds,
                                   casos.n = casos.n,
                                   theta = casos.par[1, ])) # 36 min approx
system.time(testing.p2 <- mclapply(1:4000, prob,
                                   seeds = seeds,
                                   casos.n = casos.n,
                                   theta = casos.par[2, ])) # 34 min approx

setwd("./results.RD/")

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
