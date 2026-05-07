Ctau <- function(quantile = 0.5) {
	dir <- "~/Documents/ProyectoTitulo/R/gamlss/"
	if (quantile > 0 && quantile < 1) {
		write.table(quantile, paste(dir, "tau", sep = ""),
								col.names = FALSE, row.names = FALSE)
		tau <- read.table(paste(dir, "tau", sep = ""))
		print(paste("Tau:", tau))
		return(tau[1, 1])
	}
	else
		print("Error.")
}
