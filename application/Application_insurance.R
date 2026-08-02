library(gamlss)
library(VGAM)
library(dplyr)
library(ggplot2)
library(car)

source("./aux/dBS2.R"); source("./aux/BS2.R")
source("./../models/RD.R")

insurance <- read.table("insurance.txt", header = TRUE, dec = ".")

insurance <- insurance %>%
	filter(legrep == 1) %>%
	select(-legrep)

Rd <- gamlss(amount ~ optime, sigma.fo = ~optime,
						 data = insurance,
						 family = RD, method = RS(100))
Gamma <- gamlss(amount ~ optime, sigma.fo = ~optime,
								data = insurance,
								family = GA, method = RS(100))
Wei <- gamlss(amount ~ optime, sigma.fo = ~optime,
							data = insurance,
							family = WEI3, method = RS(100))
Igamma <- gamlss(amount ~ optime, sigma.fo = ~optime,
								 data = insurance,
								 family = IG, method = RS(100))
BS <- gamlss(amount ~ optime, sigma.fo = ~optime,
						 data = insurance,
						 family = BS2, method = RS(100))

AICs <- list(Rd = Rd$aic, Gamma = Gamma$aic, Igamma = Igamma$aic,
						 Wei = Wei$aic, BS = BS$aic)
BICs <- list(Rd = Rd$sbc, Gamma = Gamma$sbc, Igamma = Igamma$sbc,
						 Wei = Wei$sbc, BS = BS$sbc)
