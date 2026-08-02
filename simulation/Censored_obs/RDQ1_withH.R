library(gamlss)
library(gamlss.cens)
library(VGAM)

RDQ1 <- function (mu.link="log", sigma.link="log", nu.link="log") {
    mstats <- checklink("mu.link", "Reparametrized.Dagum1", substitute(mu.link), c("inverse", "log", "identity"))# dummy
    dstats <- checklink("sigma.link", "Reparametrized.Dagum1", substitute(sigma.link), c("inverse", "log", "identity"))
    vstats <- checklink("nu.link", "Reparametrized.Dagum1", substitute(nu.link), c("inverse", "log", "identity"))

    structure(
          list(family = c("RDQ1", "Reparametrized.Dagum1"),
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
                    tau = 0.5
                    return(sigma*(nu*tau^(1/nu) - nu + tau^(1/nu)*(y/mu)^sigma)/(mu*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)))
                },
                 dldd = function(y, mu, sigma, nu) {
                    tau = 0.5
                    return((nu*sigma*(-log(mu) + log(y))*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1) - sigma*tau^(1/nu)*(y/mu)^sigma*(nu + 1)*log(y/mu) + tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)/(sigma*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)))
                },
                 dldv = function(y, mu, sigma, nu) {
                    tau = 0.5
                    return((nu^2*(tau^(1/nu) - 1)*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)*(-sigma*log(mu) + sigma*log(y) - log((1 - tau^(1/nu))/tau^(1/nu)) - log((-tau^(1/nu)*(y/mu)^sigma + tau^(1/nu) - 1)/(tau^(1/nu) - 1))) + nu*(tau^(1/nu) - 1)*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1) + nu*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)*log(tau) - tau^(1/nu)*(y/mu)^sigma*(nu + 1)*log(tau))/(nu^2*(tau^(1/nu) - 1)*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)))
                },
               d2ldm2 = function(y, mu, sigma, nu) {
                tau = 0.5
                return(sigma*(sigma*tau^(1/nu)*(y/mu)^sigma*(nu*tau^(1/nu) - nu + tau^(1/nu)*(y/mu)^sigma) + (tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)*(-nu*tau^(1/nu) + nu - sigma*tau^(1/nu)*(y/mu)^sigma - tau^(1/nu)*(y/mu)^sigma))/(mu^2*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)^2))
            },
              d2ldmdd = function(y, mu, sigma, nu) {
                tau = 0.5
                return((-sigma*tau^(1/nu)*(y/mu)^sigma*(nu*tau^(1/nu) - nu + tau^(1/nu)*(y/mu)^sigma)*log(y/mu) + (tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)*(nu*tau^(1/nu) - nu + sigma*tau^(1/nu)*(y/mu)^sigma*log(y/mu) + tau^(1/nu)*(y/mu)^sigma))/(mu*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)^2))
            },
              d2ldmdv = function(y, mu, sigma, nu) {
                tau = 0.5
                return(sigma*(tau^(1/nu)*((y/mu)^sigma - 1)*(nu*tau^(1/nu) - nu + tau^(1/nu)*(y/mu)^sigma)*log(tau) - (nu^2*(1 - tau^(1/nu)) + nu*tau^(1/nu)*log(tau) + tau^(1/nu)*(y/mu)^sigma*log(tau))*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1))/(mu*nu^2*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)^2))
            },
               d2ldd2 = function(y, mu, sigma, nu) {
                tau = 0.5
                return((nu*sigma^2*tau^(2/nu)*(y/mu)^sigma*log(y/mu)^2 - nu*sigma^2*tau^(1/nu)*(y/mu)^sigma*log(y/mu)^2 + sigma^2*tau^(2/nu)*(y/mu)^sigma*log(y/mu)^2 - sigma^2*tau^(1/nu)*(y/mu)^sigma*log(y/mu)^2 - tau^(2/nu)*(y/mu)^(2*sigma) + 2*tau^(2/nu)*(y/mu)^sigma - tau^(2/nu) - 2*tau^(1/nu)*(y/mu)^sigma + 2*tau^(1/nu) - 1)/(sigma^2*(tau^(2/nu)*(y/mu)^(2*sigma) - 2*tau^(2/nu)*(y/mu)^sigma + tau^(2/nu) + 2*tau^(1/nu)*(y/mu)^sigma - 2*tau^(1/nu) + 1)))
            },
              d2ldddv = function(y, mu, sigma, nu) {
                tau = 0.5
                return((-tau^(1/nu)*((y/mu)^sigma - 1)*(nu*sigma*(log(mu) - log(y))*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1) + sigma*tau^(1/nu)*(y/mu)^sigma*(nu + 1)*log(y/mu) - tau^(1/nu)*(y/mu)^sigma + tau^(1/nu) - 1)*log(tau) + (tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)*(-nu^2*sigma*(tau^(1/nu)*(y/mu)^sigma*log(y/mu) + (log(mu) - log(y))*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)) + nu*sigma*tau^(1/nu)*((y/mu)^sigma - 1)*(log(mu) - log(y))*log(tau) + tau^(1/nu)*(sigma*(y/mu)^sigma*(nu + 1)*log(y/mu) - (y/mu)^sigma + 1)*log(tau)))/(nu^2*sigma*(tau^(1/nu)*(y/mu)^sigma - tau^(1/nu) + 1)^2))
            },
            d2ldv2 = function(y, mu, sigma, nu) {
                tau = 0.5
                return((-nu^2*tau^(4/nu)*(y/mu)^(2*sigma) + 2*nu^2*tau^(4/nu)*(y/mu)^sigma - nu^2*tau^(4/nu) + 2*nu^2*tau^(3/nu)*(y/mu)^(2*sigma) - 6*nu^2*tau^(3/nu)*(y/mu)^sigma + 4*nu^2*tau^(3/nu) - nu^2*tau^(2/nu)*(y/mu)^(2*sigma) + 6*nu^2*tau^(2/nu)*(y/mu)^sigma - 6*nu^2*tau^(2/nu) - 2*nu^2*tau^(1/nu)*(y/mu)^sigma + 4*nu^2*tau^(1/nu) - nu^2 + 2*nu*tau^(3/nu)*(y/mu)^(2*sigma)*log(tau) - nu*tau^(3/nu)*(y/mu)^sigma*log(tau)^2 - 2*nu*tau^(3/nu)*(y/mu)^sigma*log(tau) + nu*tau^(3/nu)*log(tau)^2 - 2*nu*tau^(2/nu)*(y/mu)^(2*sigma)*log(tau) + 2*nu*tau^(2/nu)*(y/mu)^sigma*log(tau)^2 + 4*nu*tau^(2/nu)*(y/mu)^sigma*log(tau) - 2*nu*tau^(2/nu)*log(tau)^2 - nu*tau^(1/nu)*(y/mu)^sigma*log(tau)^2 - 2*nu*tau^(1/nu)*(y/mu)^sigma*log(tau) + nu*tau^(1/nu)*log(tau)^2 - tau^(3/nu)*(y/mu)^(2*sigma)*log(tau)^2 + tau^(3/nu)*(y/mu)^sigma*log(tau)^2 - tau^(1/nu)*(y/mu)^sigma*log(tau)^2)/(nu^4*(tau^(4/nu)*(y/mu)^(2*sigma) - 2*tau^(4/nu)*(y/mu)^sigma + tau^(4/nu) - 2*tau^(3/nu)*(y/mu)^(2*sigma) + 6*tau^(3/nu)*(y/mu)^sigma - 4*tau^(3/nu) + tau^(2/nu)*(y/mu)^(2*sigma) - 6*tau^(2/nu)*(y/mu)^sigma + 6*tau^(2/nu) + 2*tau^(1/nu)*(y/mu)^sigma - 4*tau^(1/nu) + 1)))
            },
          G.dev.incr  = function(y,mu,sigma,nu,...) -2*dRDQ1(y, mu, sigma, nu, log = TRUE), 
                rqres = expression(rqres(pfun="pRDQ1", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
           mu.initial = expression(mu <- rep(median(y),length(y))), 
        sigma.initial = expression(sigma <- rep(1, length(y))), 
           nu.initial = expression(nu <- rep(1, length(y))),
             mu.valid = function(mu) all(mu > 0), 
          sigma.valid = function(sigma)  all(sigma > 0),
             nu.valid = function(nu) all(nu > 0),
              y.valid = function(y)  all(y > 0),
                 mean = function(mu, sigma, nu) {
                    a = sigma
                    b = mu * (0.5^(-1/nu)-1)^(1/sigma)
                    p = nu
                    return(ifelse(a > 1, p * b * beta(p+1/a, 1-1/a), Inf))
                  },
             variance = function(mu, sigma, nu) {
                a = sigma
                b = mu * (0.5^(-1/nu)-1)^(1/sigma)
                p = nu
                return(ifelse(a > 2, b^2 * (beta(1-2/a, p+2/a) - beta(1-1/a, p+1/a)^2), Inf))
              }
          ),
                class = c("gamlss.family","family"))
}

dRDQ1 <- function(x, mu = 1, sigma = 1, nu = 1.5, log = FALSE) {
  tau <- 0.5
  mu2 <- mu * (tau^(-1 / nu) - 1)^(1 / sigma)
  ddagum(x = x, scale = mu2, shape1.a = sigma, shape2.p = nu, log = log)
}

pRDQ1 <- function(q, mu = 1, sigma = 1, nu = 1.5, lower.tail = TRUE, log.p = FALSE) {
  tau <- 0.5
  mu2 <- mu * (tau^(-1 / nu) - 1)^(1 / sigma)
  pdagum(q = q, scale = mu2, shape1.a = sigma, shape2.p = nu, lower.tail = lower.tail, log.p = log.p)
}

qRDQ1 <- function(p, mu = 1, sigma = 1, nu = 1.5, lower.tail = TRUE, log.p = FALSE) {
  tau <- 0.5
  mu2 <- mu * (tau^(-1 / nu) - 1)^(1 / sigma)
  qdagum(p = p, scale = mu2, shape1.a = sigma, shape2.p = nu, lower.tail = lower.tail, log.p = log.p)
}

rRDQ1 <- function(n, mu = 1, sigma = 1, nu = 1.5) {
  tau <- 0.5
  mu2 <- mu * (tau^(-1 / nu) - 1)^(1 / sigma)
  rdagum(n = n, scale = mu2, shape1.a = sigma, shape2.p = nu)
}

# Log Verosimilitud

library(pracma)

LogLS <- function(theta, DM, data) {
  mu_ <- theta[1:3]; sigma_ <- theta[4:5]; nu_ <- theta[6]

  X1 <- DM[, 1:3]; X2 <- DM[, 4:5]; X3 <- DM[, 6]

  mu <- as.vector(exp(X1 %*% mu_))        # log-link
  sigma <- as.vector(exp(X2 %*% sigma_))  # log-link
  nu <- as.vector(exp(X3 * nu_))          # log-link-desplazado (nu > 1)

  T <- data[, 1]          # T: Tiempo observado
  I <- data[, 2]          # I: indicadora falla/censura
  indx_C <- which(I == 0) # Indice elementos censurados

  T_C <- T[indx_C]
  mu_C <- mu[indx_C]; sigma_C <- sigma[indx_C]; nu_C <- nu[indx_C]

  T_F <- T[-indx_C]
  mu_F <- mu[-indx_C]; sigma_F <- sigma[-indx_C]; nu_F <- nu[-indx_C]

  return(sum(dRDQ1(T_F, mu = mu_F, sigma = sigma_F, nu = nu_F, log = TRUE)) +
         sum(pRDQ1(T_C, mu = mu_C, sigma = sigma_C, nu = nu_C, lower.tail = FALSE, log.p = TRUE)))
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
  y <- rRDQ1(n, mu.true, sigma.true, nu.true)  # T: tiempos de falla
  lmbds <- apply(cbind(mu.true, sigma.true, nu.true), 1, find_delta,
                 model = dRDQ1, censorship = censorship, search_interval = c(0, 100))
  C <- rexp(n = n, lmbds)
  
  # GAMLSS
  aux <- gamlss(Surv(time = pmin(y, C), event = y < C) ~ X1[, -1, drop=FALSE],
                sigma.fo = ~X2[, -1, drop = FALSE],
                family = cens(RDQ1), method = RS(1000),
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
    print(errors)
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
                                     censorship = Cens)) # 63 min approx cens = 0.1
  system.time(testing.p2 <- mclapply(1:4000, prob,
                                     seeds = seeds,
                                     casos.n = casos.n,
                                     theta = casos.par[2, ],
                                     censorship = Cens)) # 77 min approx cens = 0.1
  
  setwd(paste("./RDQ1/RDQ1_", Cens_level, sep = ""))
  
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
