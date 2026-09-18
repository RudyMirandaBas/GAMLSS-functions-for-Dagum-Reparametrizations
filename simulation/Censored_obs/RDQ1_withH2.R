## ============================================================
## Uso: Rscript MC_RDQ1_censored_simulation.R <tau> <censorship>
## Ejemplo: Rscript MC_RDQ1_censored_simulation.R 0.5 0.25
## ============================================================

## --- Leer argumentos desde la linea de comandos ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Debe indicar tau y censorship como argumentos.\n",
       "Uso: Rscript MC_RDQ1_censored_simulation.R <tau> <censorship>\n",
       "Ejemplo: Rscript MC_RDQ1_censored_simulation.R 0.5 0.25", call. = FALSE)
}

tau <- as.numeric(args[1])
censorship <- as.numeric(args[2])

if (is.na(tau) || tau <= 0 || tau >= 1) {
  stop("tau debe ser un numero entre 0 y 1 (ej: 0.1, 0.25, 0.5, 0.75, 0.9). Valor recibido: ", args[1], call. = FALSE)
}

if (is.na(censorship) || censorship <= 0 || censorship >= 1) {
  stop("censorship debe ser un numero entre 0 y 1 (ej: 0.1, 0.25, 0.5). Valor recibido: ", args[2], call. = FALSE)
}

cat("Corriendo simulacion con tau =", tau, " y censorship =", censorship, "\n")

library(gamlss)
library(gamlss.cens)
library(VGAM)

source("./../../models/RDQ1.R")

# Log Verosimilitud

library(pracma)

LogLS <- function(theta, DM, data, tau) {
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

  return(sum(dRDQ1(T_F, mu = mu_F, sigma = sigma_F, nu = nu_F, tau = tau, log = TRUE)) +
         sum(pRDQ1(T_C, mu = mu_C, sigma = sigma_C, nu = nu_C, tau = tau, lower.tail = FALSE, log.p = TRUE)))
}

#######################################
############# Monte Carlo #############
#######################################

library(parallel)
source("./find_delta.R")

MCllikelihood_estimation <- function(n = 1000, mu_ = c(1, .5, .2), sigma_ = c(.5, .4), nu_ = c(-.2), tau = 0.5, censorship = 0.1) {

  z1 <- as.numeric(scale(runif(n)))
  z2 <- as.numeric(scale(rnorm(n)))

  X1 <- model.matrix(~z1+z2)                   #matriz de diseño de mu
  X2 <- model.matrix(~z1)                      #matriz de diseño de sigma
  X3 <- matrix(1, nrow = n)                    #matriz de diseño de nu
  colnames(X3) <- "(Intercept)"

  mu.true <- as.vector(exp(X1 %*% mu_))        # log-link
  sigma.true <- as.vector(exp(X2 %*% sigma_))  # log-link
  nu.true <- as.vector(exp(X3 %*% nu_))        # log-link-desplazado (nu > 1)
  y <- rRDQ1(n, mu.true, sigma.true, nu.true, tau)  # T: tiempos de falla
  lmbds <- apply(cbind(mu.true, sigma.true, nu.true), 1, find_delta,
                 model = dRDQ1, censorship = censorship, search_interval = c(0, 100))
  C <- rexp(n = n, lmbds)

  # GAMLSS
  aux <- gamlss(Surv(time = pmin(y, C), event = y < C) ~ X1[, -1, drop=FALSE],
                sigma.fo = ~X2[, -1, drop = FALSE],
                tau.fix = TRUE, tau.start = tau,
                family = cens(RDQ1), method = RS(1000),
                control = gamlss.control(trace = FALSE))

  res.gamlss <- c(
    coef(aux, what = "mu"),
    coef(aux, what = "sigma"),
    coef(aux, what = "nu")
  )

  # Solo estimación puntual; SE de GAMLSS no se calcula
  res <- cbind(
    c(mu_, sigma_, nu_),   # REAL
    res.gamlss,            # estimación GAMLSS
    NA                     # SE GAMLSS
  )

  H <- hessian(LogLS, x0 = res[, 2], DM = cbind(X1, X2, X3),
               data = Surv(time = pmin(y, C), event = y < C),
               tau = tau)

  VCOV <- sqrt(diag(solve(-H)))

  res <- cbind(res, VCOV)

  colnames(res) <- c("true",
                     "est gamlss", "se gamlss", "se hessian")
  rownames(res) <- c(paste("beta1", 1:ncol(X1), sep = ""),
                     paste("beta2", 1:ncol(X2), sep = ""),
                     paste("beta3", 1:ncol(X3), sep = ""))

  return(list(Results = res,
              LL = c(gamlss = logLik(aux)),
              Converged = aux$converged,
              Iterations = aux$iter))       # <-- iteraciones de gamlss
}

# to use in mclapply
RNGkind("L'Ecuyer-CMRG")

set.seed(123)
seeds <- list(.Random.seed)
for (i in 2:4000) {
  seeds[[i]] <- nextRNGStream(seeds[[i - 1]])
}

prob <- function(i, theta, seeds, casos.n, tau, censorship) {
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
  elapsed_secs <- NA          # <-- tiempo del intento exitoso (segundos)

  while (temp) {
    t_start <- Sys.time()     # <-- inicio de este intento

    testing <- try(MCllikelihood_estimation(n = n,
                                            mu_ = theta[1:3],
                                            sigma_ = theta[4:5],
                                            nu_ = theta[6],
                                            tau = tau,
                                            censorship = censorship),
                   silent = TRUE)

    t_end <- Sys.time()       # <-- fin de este intento

    temp <- grepl("Error", testing)[1]
    errors <- errors + 1

    if (!temp) {
      # este intento fue exitoso -> guardamos su tiempo
      elapsed_secs <- as.numeric(difftime(t_end, t_start, units = "secs"))
    }
  }

  testing[["Errors"]] <- errors
  testing[["Time"]] <- elapsed_secs   # <-- tiempo (seg) del intento exitoso
  # testing[["Iterations"]] ya viene incluido desde MCllikelihood_estimation()

  cat("%:", round(i / 4000 * 100, 3), " (", i, "/", 4000, ")  n:", n,
      "  time:", round(elapsed_secs, 3), "s  gamlss iters:", testing[["Iterations"]],
      "\n", sep = "")
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
                                   theta = casos.par[1, ],
                                   tau = tau,
                                   censorship = censorship)) # 63 min approx cens = 0.1
#system.time(testing.p2 <- mclapply(1:4000, prob,
#                                   seeds = seeds,
#                                   casos.n = casos.n,
#                                   theta = casos.par[2, ],
#                                   tau = tau,
#                                   censorship = censorship)) # 77 min approx cens = 0.1

## --- Carpeta de resultados: RDQ1 / RDQ1_<censura> / tau_<tau> ---
cens_pct <- censorship * 100
cens_str <- if (isTRUE(all.equal(cens_pct, round(cens_pct)))) {
  as.character(as.integer(round(cens_pct)))                    # 0.1 -> "10", 0.25 -> "25"
} else {
  gsub("\\.", "_", format(cens_pct, trim = TRUE))               # 0.125 -> "12_5"
}

tau_str <- gsub("\\.", "_", format(tau, trim = TRUE))           # 0.5 -> "0_5"

results_dir <- file.path(".", "RDQ1", paste0("RDQ1_", cens_str), paste0("tau_", tau_str))

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
setwd(results_dir)

reps <- 1000
for (i in 1:length(casos.n)) {

  ############ set 1

  temp_p1 <- testing.p1[(1 + reps * (i - 1)):(reps * i)]

  # LL
  LL_p1 <- sapply(temp_p1, "[[", "LL")
	name_ll_p1 <- paste("LL_p1_", casos.n[i], '.csv', sep = '')
  write.csv(LL_p1, file = name_ll_p1, row.names = FALSE)

  # pars
  pars_p1 <- t(sapply(temp_p1, "[[", "Results"))
	name_pars_p1 <- paste("pars_p1_", casos.n[i], '.csv', sep = '')
  write.csv(pars_p1, file = name_pars_p1, row.names = FALSE)

  # Converged
  conv_p1 <- sapply(temp_p1, "[[", "Converged")
	name_conv_p1 <- paste("noconv_p1_", casos.n[i], '.csv', sep = '')
  write.csv(sum(!conv_p1), file = name_conv_p1, row.names = FALSE)

  # Errors
  errors_p1 <- t(sapply(temp_p1, "[[", "Errors"))
	name_errors_p1 <- paste("errors_p1_", casos.n[i], '.csv', sep = '')
  write.csv(sum(errors_p1), file = name_errors_p1, row.names = FALSE)

  # Time (tiempo en segundos del intento exitoso)
  time_p1 <- sapply(temp_p1, "[[", "Time")
	name_time_p1 <- paste("time_p1_", casos.n[i], '.csv', sep = '')
  write.csv(time_p1, file = name_time_p1, row.names = FALSE)

  # Iterations (numero de iteraciones de gamlss)
  giter_p1 <- sapply(temp_p1, "[[", "Iterations")
	name_giter_p1 <- paste("giter_p1_", casos.n[i], '.csv', sep = '')
  write.csv(giter_p1, file = name_giter_p1, row.names = FALSE)

#  ############ set 2
#
#  temp_p2 <- testing.p2[(1 + reps * (i - 1)):(reps * i)]
#
#  # LL
#  LL_p2 <- sapply(temp_p2, "[[", "LL")
#	name_ll_p2 <- paste("LL_p2_", casos.n[i], '.csv', sep = '')
#  write.csv(LL_p2, file = name_ll_p2, row.names = FALSE)
#
#  # pars
#  pars_p2 <- t(sapply(temp_p2, "[[", "Results"))
#	name_pars_p2 <- paste("pars_p2_", casos.n[i], '.csv', sep = '')
#  write.csv(pars_p2, file = name_pars_p2, row.names = FALSE)
#
#  # Converged
#  conv_p2 <- sapply(temp_p2, "[[", "Converged")
#	name_conv_p2 <- paste("noconv_p2_", casos.n[i], '.csv', sep = '')
#  write.csv(sum(!conv_p2), file = name_conv_p2, row.names = FALSE)
#
#  # Errors
#  errors_p2 <- t(sapply(temp_p2, "[[", "Errors"))
#	name_errors_p2 <- paste("errors_p2_", casos.n[i], '.csv', sep = '')
#  write.csv(sum(errors_p2), file = name_errors_p2, row.names = FALSE)
#
#  # Time (tiempo en segundos del intento exitoso)
#  time_p2 <- sapply(temp_p2, "[[", "Time")
#	name_time_p2 <- paste("time_p2_", casos.n[i], '.csv', sep = '')
#  write.csv(time_p2, file = name_time_p2, row.names = FALSE)
#
#  # Iterations (numero de iteraciones de gamlss)
#  giter_p2 <- sapply(temp_p2, "[[", "Iterations")
#	name_giter_p2 <- paste("giter_p2_", casos.n[i], '.csv', sep = '')
#  write.csv(giter_p2, file = name_giter_p2, row.names = FALSE)
}

cat("Simulacion terminada. Resultados guardados en:", normalizePath(getwd()), "\n")
