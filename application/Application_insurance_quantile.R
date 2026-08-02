library(gamlss)
library(dplyr)
library(VGAM)
library(car)
library(alr4); attach(landrent)

source("./../models/RDQ1.R")
source("./../models/RDQ2.R")

plot(X2, Y / X1, ylab = "ratio(y)", xlab = "density of dairy cows (x)", pch = 19)

snacks_RDQ1_gl <- suppressWarnings(try(gamlss(formula = I(Y / X1) ~ X2,
                                              sigma.fo = ~ 1,
                                              tau.fix = TRUE,
                                              tau.start = .4,
			                                        family = RDQ1,
														                  method = RS(1000),
                                              data = landrent)))
summary(snacks_RDQ1_gl)

snacks_RDQ2_gl <- suppressWarnings(try(gamlss(formula = I(Y / X1) ~ X2,
                                              sigma.fo = ~ 1,
                                              tau.fix = TRUE,
                                              tau.start = .4,
			                                        family = RDQ2,
														                  method = RS(1000),
                                              data = landrent)))
summary(snacks_RDQ2_gl)

logLik(snacks_RDQ1_gl)
logLik(snacks_RDQ2_gl)

ord <- order(X2)

pdf("RDQ1and2.pdf")
par(mai = c(0.85,0.99,0.05,0.05))
plot(X2, Y / X1, ylab = "ratio(y)", xlab = "density of dairy cows (x)", pch = 19,
     cex.axis = 1.6, cex.lab = 1.8)
lines(X2[ord], fitted(snacks_RDQ1_gl, "mu")[ord], col = "blue", lwd = 2)
lines(X2[ord], fitted(snacks_RDQ2_gl, "mu")[ord], col = "red", lwd = 2)
legend("topright",
       legend = c("RDQ1", "RDQ2"),
       col = c("blue", "red"), lty = 1, lwd = 2, bty = "n", cex = 1.5)
dev.off()

# Residuos

pdf("Residuals_RDQ1.pdf", width = 5, height = 5, pointsize = 7)
par(mai = c(0.85,0.99,0.05,0.05))
par(cex.axis = 1.7, cex.lab = 2)
qqPlot(snacks_RDQ1_gl$residuals, envelope = list(style = "lines"),
       col.lines = "black", grid = FALSE, id = FALSE, pch = 21, cex = 1,
       ylab = "Empiric Quantiles")
dev.off()

pdf("Residuals_RDQ2.pdf", width = 5, height = 5, pointsize = 7)
par(mai = c(0.85,0.99,0.05,0.05))
par(cex.axis = 1.7, cex.lab = 2)
qqPlot(snacks_RDQ2_gl$residuals, envelope = list(style = "lines"),
       col.lines = "black", grid = FALSE, id = FALSE, pch = 21, cex = 1,
       ylab = "Empiric Quantiles")
dev.off()
