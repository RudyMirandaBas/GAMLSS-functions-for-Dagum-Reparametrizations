model <- c("RD", "RDQ1", "RDQ2")
dir <- paste("~/Documents/ProyectoTitulo/R/MC/Censored/RDQ2/HESSIAN/RDQ2_50/")
             #model[3], sep = "")
setwd(dir)

casos.n <- c(50, 100, 200, 500)
casos.par <- matrix(c(1, .5, .2, .5, .4, -.2, -1, -.5, -.2, .25, .7, .3),
                    nrow = 2, byrow = TRUE)

### LogLikelihoods

likelihoods <- t(t(sapply(lapply(paste("LL_p",
                                       rep(1:2, each = 4), "_",
                                       rep(casos.n, 2), ".csv",
                                       sep = ""),
                                 read.csv), colMeans)))
rownames(likelihoods) <- paste("B", rep(1:2, each = 4), "_", casos.n, sep = "")

### Errors

matrix(unlist(lapply(paste("errors_p", rep(1:2, each = 4), "_", rep(casos.n, 2), ".csv", sep = ""), read.csv)),
       nr = 4, dimnames = list(casos.n, c("p1", "p2")))

### MSE
MSE <- lapply(lapply(paste("pars_p", rep(1:2, each = 4), "_", rep(casos.n, 2), ".csv", sep = ""), read.csv), function(x) {
       vec <- sqrt(colSums((x[1:6] - x[7:12])^2) / dim(x)[1])
       names(vec) <- c(paste("beta1", 1:3, sep = ""),
                       paste("beta2", 1:2, sep = ""),
                       paste("beta3", 1:1, sep = ""))
       return(vec)
  }
)
names(MSE) <- paste("B", rep(1:2, each = 4), "_", casos.n, sep = "")

### Coverage Probability GAMLSS
CP <- lapply(lapply(paste("pars_p", rep(1:2, each = 4), "_", rep(casos.n, 2), ".csv", sep = ""), read.csv), function(x) {
       b <- x[1:6]                                # Valor Real de Beta
       se <- x[13:18]                             # se gamlss
       b_hat <- x[7:12]                           # Valor estimado de Beta gamlss
       lb <- b_hat - qnorm(1 - .05 / 2) * se      # Lower Bound CI
       ub <- b_hat + qnorm(1 - .05 / 2) * se      # Upper Bound CI
       vec <- colMeans(ifelse(b > lb & b < ub, 1, 0))

       names(vec) <- c(paste("beta1", 1:3, sep = ""),
                       paste("beta2", 1:2, sep = ""),
                       paste("beta3", 1:1, sep = ""))
       return(vec)
  }
)
names(CP) <- paste("B", rep(1:2, each = 4), "_", casos.n, sep = "")

### Coverage Probability HESSIAN
CPH <- lapply(lapply(paste("pars_p", rep(1:2, each = 4), "_", rep(casos.n, 2), ".csv", sep = ""), read.csv), function(x) {
       b <- x[1:6]                                # Valor Real de Beta
       se <- x[19:24]                             # se gamlss
       b_hat <- x[7:12]                           # Valor estimado de Beta gamlss
       lb <- b_hat - qnorm(1 - .05 / 2) * se      # Lower Bound CI
       ub <- b_hat + qnorm(1 - .05 / 2) * se      # Upper Bound CI
       vec <- colMeans(ifelse(b > lb & b < ub, 1, 0), na.rm = TRUE)

       names(vec) <- c(paste("beta1", 1:3, sep = ""),
                       paste("beta2", 1:2, sep = ""),
                       paste("beta3", 1:1, sep = ""))
       return(vec)
  }
)
names(CPH) <- paste("B", rep(1:2, each = 4), "_", casos.n, sep = "")

### Betas

matrix_names <- list(c(paste("beta1", 1:3, sep = ""),
                       paste("beta2", 1:2, sep = ""),
                       paste("beta3", 1:1, sep = "")),
                     c("true", "est_gamlss", "se_gamlss", "se_hessian"))

Coefficients <- lapply(lapply(paste("pars_p", rep(1:2, each = 4), "_", rep(casos.n, 2), ".csv", sep = ""), read.csv),
       function(x) matrix(colMeans(x, na.rm = TRUE), nr = 6, dimnames = matrix_names))

names(Coefficients) <- paste("B", rep(1:2, each = 4), "_", casos.n, sep = "")

for (i in 1:length(Coefficients)) {
  x <- Coefficients[[i]]
  Coefficients[[i]] <- cbind(x, Bias = x[, 2] - x[, 1], RMSE = MSE[[i]], CP = CP[[i]], CPH = CPH[[i]])
  Coefficients[[i]] <- Coefficients[[i]][, c(5, 3, 6, 8)]
}

print(Coefficients)
