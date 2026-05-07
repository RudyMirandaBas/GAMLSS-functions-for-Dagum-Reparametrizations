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
                b = exp(log(mu) - log(sigma) - lbeta(sigma+1/nu, 1-1/nu))
                p = sigma
                return(ifelse(a > 2, b^2 * (beta(1-2/a, p+2/a) - beta(1-1/a, p+1/a)^2), Inf))
              }
          ),
                class = c("gamlss.family","family"))
}

dRD <- function(x, mu = 1, sigma = 1, nu = 1.5, log = FALSE) {
  mu2 <- exp(log(mu) - log(sigma) - lbeta(sigma+1/nu, 1-1/nu))
  fx <- ddagum(x = x, scale = mu2, shape1.a = nu, shape2.p = sigma, log = log)
  ifelse(x > 0, fx, 0)
}

pRD <- function(q, mu = 1, sigma = 1, nu = 1.5, lower.tail = TRUE, log.p = FALSE) {
  mu2 <- exp(log(mu) - log(sigma) - lbeta(sigma+1/nu, 1-1/nu))
  pdagum(q = q, scale = mu2, shape1.a = nu, shape2.p = sigma, lower.tail = lower.tail, log.p = log.p)
}

qRD <- function(p, mu = 1, sigma = 1, nu = 1.5, lower.tail = TRUE, log.p = FALSE) {
  mu2 <- exp(log(mu) - log(sigma) - lbeta(sigma+1/nu, 1-1/nu))
  qdagum(p = p, scale = mu2, shape1.a = nu, shape2.p = sigma, lower.tail = lower.tail, log.p = log.p)
}

rRD <- function(n, mu = 1, sigma = 1, nu = 1.5) {
  mu2 <- exp(log(mu) - log(sigma) - lbeta(sigma+1/nu, 1-1/nu))
  rdagum(n = n, scale = mu2, shape1.a = nu, shape2.p = sigma)
}

#dat <- rRD(1000, 5, 4, 1.8)
#gamlss(dat~1, family=RD)
#gamlss(dat~1, family=WEI)
#
#l <- function(theta, x) {
#  mu=theta[1]; sigma=theta[2]; nu=theta[3]
#  # Define B
#  B = beta(sigma+1/nu, 1-1/nu)
#
#  # Define L
#  (1+sigma*nu)*log(sigma) + log(nu) + (sigma*nu-1)*log(x) - sigma*nu*log(mu) + sigma*nu*log(B) - (sigma+1)*log(1+(x*sigma*B/mu)^nu)
#}
#
#optim(par = c(5,1,1.5), fn = function(x) -sum(dRD(x = dat, mu = x[1], sigma = x[2], nu = x[3], log = TRUE)), control=list(maxit=10000))
