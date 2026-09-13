## ============================================================
## Uso: Rscript MC_RDQ2_simulation.R <tau>
## Ejemplo: Rscript MC_RDQ2_simulation.R 0.5
## ============================================================

## --- Leer argumento tau desde la linea de comandos ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Debe indicar el valor de tau como argumento.\n",
       "Uso: Rscript MC_RDQ2_simulation.R <tau>\n",
       "Ejemplo: Rscript MC_RDQ2_simulation.R 0.5", call. = FALSE)
}

tau <- as.numeric(args[1])

if (is.na(tau) || tau <= 0 || tau >= 1) {
  stop("tau debe ser un numero entre 0 y 1 (ej: 0.5). Valor recibido: ", args[1], call. = FALSE)
}

cat("Corriendo simulacion con tau =", tau, "\n")

library(gamlss)
library(VGAM)

source("./../../models/RDQ2.R")

#######################################
############# Monte Carlo #############
#######################################

library(parallel)

MCllikelihood_estimation <- function(n = 1000, mu_ = c(1, .5, .2), sigma_ = c(.5, .4), nu_ = c(-.2), tau = 0.5) {

  z1 <- as.numeric(scale(runif(n)))
  z2 <- as.numeric(scale(rnorm(n)))

  X1 <- model.matrix(~z1+z2)                   #matriz de diseño de mu
  X2 <- model.matrix(~z1)                      #matriz de diseño de sigma
  X3 <- matrix(1, nrow = n)                    #matriz de diseño de nu
  colnames(X3) <- "(Intercept)"

  mu.true <- as.vector(exp(X1 %*% mu_))        # log-link
  sigma.true <- as.vector(exp(X2 %*% sigma_))  # log-link
  nu.true <- as.vector(exp(X3 %*% nu_))        # log-link
  y <- rRDQ2(n, mu.true, sigma.true, nu.true, tau)  # simulando valores

  # GAMLSS
  aux <- gamlss(y ~ X1[, -1, drop=FALSE],
                sigma.fo = ~X2[, -1, drop = FALSE],
                tau.fix = TRUE, tau.start = tau,
                family = RDQ2, method = RS(1000),
                control = gamlss.control(trace = FALSE))
  capture.output(res.gamlss <- summary(aux)[, 1:2], file = nullfile())

  res <- cbind(c(mu_, sigma_, nu_),                              # REAL
               res.gamlss[-7, ])                                 # GAMLSS

  colnames(res) <- c("true",
                     "est gamlss", "se gamlss")
  rownames(res) <- c(paste("beta1", 1:ncol(X1), sep = ""),
                     paste("beta2", 1:ncol(X2), sep = ""),
                     paste("beta3", 1:ncol(X3), sep = ""))
  return(list(Results = res,
              LL = c(gamlss = logLik(aux)),
              Converged = aux$converged,
              Iterations = aux$iter))
}

# to use in mclapply
RNGkind("L'Ecuyer-CMRG")

set.seed(123)
seeds <- list(.Random.seed)
for (i in 2:4000) {
  seeds[[i]] <- nextRNGStream(seeds[[i - 1]])
}

prob <- function(i, theta, seeds, casos.n, tau) {
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
  elapsed_secs <- NA

  while (temp) {
    t_start <- Sys.time()

    testing <- try(MCllikelihood_estimation(n = n,
                                            mu_ = theta[1:3],
                                            sigma_ = theta[4:5],
                                            nu_ = theta[6],
                                            tau = tau),
                   silent = TRUE)

    t_end <- Sys.time()

		temp <- grepl("Error", testing)[1]
    errors <- errors + 1

    if (!temp) {
      # este intento fue exitoso -> guardamos su tiempo
      elapsed_secs <- as.numeric(difftime(t_end, t_start, units = "secs"))
    }
  }

  testing[["Errors"]] <- errors
  testing[["Time"]] <- elapsed_secs

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
                                   tau = tau)) # 8 min approx
system.time(testing.p2 <- mclapply(1:4000, prob,
                                   seeds = seeds,
                                   casos.n = casos.n,
                                   theta = casos.par[2, ],
                                   tau = tau)) # 7 min approx

## --- Carpeta de resultados con el cuantil (tau) en la ruta ---
tau_str <- gsub("\\.", "_", format(tau, trim = TRUE))     # 0.5 -> "0_5", para evitar puntos en la ruta
results_dir <- paste0("./results.RDQ2/tau_", tau_str, "/")

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

  ############ set 2

  temp_p2 <- testing.p2[(1 + reps * (i - 1)):(reps * i)]

  # LL
  LL_p2 <- sapply(temp_p2, "[[", "LL")
	name_ll_p2 <- paste("LL_p2_", casos.n[i], '.csv', sep = '')
  write.csv(LL_p2, file = name_ll_p2, row.names = FALSE)

  # pars
  pars_p2 <- t(sapply(temp_p2, "[[", "Results"))
	name_pars_p2 <- paste("pars_p2_", casos.n[i], '.csv', sep = '')
  write.csv(pars_p2, file = name_pars_p2, row.names = FALSE)

  # Converged
  conv_p2 <- sapply(temp_p2, "[[", "Converged")
	name_conv_p2 <- paste("noconv_p2_", casos.n[i], '.csv', sep = '')
  write.csv(sum(!conv_p2), file = name_conv_p2, row.names = FALSE)

  # Errors
  errors_p2 <- t(sapply(temp_p2, "[[", "Errors"))
	name_errors_p2 <- paste("errors_p2_", casos.n[i], '.csv', sep = '')
  write.csv(sum(errors_p2), file = name_errors_p2, row.names = FALSE)

  # Time (tiempo en segundos del intento exitoso)
  time_p2 <- sapply(temp_p2, "[[", "Time")
	name_time_p2 <- paste("time_p2_", casos.n[i], '.csv', sep = '')
  write.csv(time_p2, file = name_time_p2, row.names = FALSE)

  # Iterations (numero de iteraciones de gamlss)
  giter_p2 <- sapply(temp_p2, "[[", "Iterations")
	name_giter_p2 <- paste("giter_p2_", casos.n[i], '.csv', sep = '')
  write.csv(giter_p2, file = name_giter_p2, row.names = FALSE)
}

cat("Simulacion terminada. Resultados guardados en:", normalizePath(getwd()), "\n")
