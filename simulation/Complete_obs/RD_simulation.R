## ============================================================
## Uso: Rscript MC_RD_simulation.R <nu>
## Ejemplo: Rscript MC_RD_simulation.R 1.05
##
## nu es el valor VERDADERO (escala natural, nu > 1) usado para
## generar los datos. Internamente se convierte al coeficiente en
## la escala del link "logshiftto1" (linkfun(nu) = log(nu - 1))
## para poder usarlo como intercepto en X3 %*% nu_.
## ============================================================
 
## --- Leer argumento nu desde la linea de comandos ---
args <- commandArgs(trailingOnly = TRUE)
 
if (length(args) < 1) {
  stop("Debe indicar el valor de nu como argumento.\n",
       "Uso: Rscript MC_RD_simulation.R <nu>\n",
       "Ejemplo: Rscript MC_RD_simulation.R 1.05", call. = FALSE)
}
 
nu <- as.numeric(args[1])
 
if (is.na(nu) || nu <= 1) {
  stop("nu debe ser un numero mayor que 1 (ej: 1.05, 1.1, 1.25, 1.5). Valor recibido: ", args[1], call. = FALSE)
}
 
cat("Corriendo simulacion con nu =", nu, "\n")
 
library(gamlss)
library(VGAM)
 
source("./../../models/RD.R")
 
#######################################
############# Monte Carlo #############
#######################################
 
library(parallel)
ncores <- max(1, detectCores() - 1)
 
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
  elapsed_secs <- NA          # <-- tiempo del intento exitoso (segundos)
 
  while (temp) {
    t_start <- Sys.time()     # <-- inicio de este intento
 
    testing <- try(MCllikelihood_estimation(n = n,
                                            mu_ = theta[1:3],
                                            sigma_ = theta[4:5],
                                            nu_ = theta[6]),
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
 
## nu_ (coeficiente en la escala del link "logshiftto1") se calcula
## a partir del nu "verdadero" (escala natural, > 1) recibido por linea
## de comandos, usando el linkfun definido en RD.R: log(nu - 1)
nu_link <- logshiftto1()
nu_coef <- nu_link$linkfun(nu)
 
casos.n <- c(50, 100, 200, 500)
casos.par <- matrix(c(1, .5, .2, .5, .4, nu_coef,
                      -1, -.5, -.2, .25, .7, nu_coef),
                    nrow = 2, byrow = TRUE)
 
# mclapply solo funciona en Linux y MacOS,
# no en windows. Si corre esto en windows hace un lapply solamente.
# Forma mas facil de aplicar calculo paralelo.
system.time(testing.p1 <- mclapply(1:4000, prob,
                                   seeds = seeds,
                                   casos.n = casos.n,
                                   theta = casos.par[1, ], mc.cores = ncores)) # 36 min approx
system.time(testing.p2 <- mclapply(1:4000, prob,
                                   seeds = seeds,
                                   casos.n = casos.n,
                                   theta = casos.par[2, ], mc.cores = ncores)) # 34 min approx
 
## --- Carpeta de resultados con el valor de nu en la ruta ---
nu_str <- gsub("\\.", "_", format(nu, trim = TRUE))     # 1.05 -> "1_05", para evitar puntos en la ruta
results_dir <- paste0("./results.RD/nu_", nu_str, "/")
 
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
