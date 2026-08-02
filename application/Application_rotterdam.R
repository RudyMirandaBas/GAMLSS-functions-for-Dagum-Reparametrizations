library(gamlss)
library(gamlss.cens)
library(VGAM)
library(dplyr)
library(tidyr)
library(survival)
library(ggplot2)
library(survminer)
library(car)

source("./../models/RD.R")
source("./aux/dBS2.R"); source("./aux/BS2.R")

# ---
q_age <- quantile(rotterdam$age, probs = seq(0, 1, .5), na.rm = TRUE)
rotterdam$age_q <- cut(rotterdam$age,
                       breaks = q_age,
                       include.lowest = TRUE,
                       labels = c("Q1", "Q2"))
rotterdam <- rotterdam %>% mutate(size = as.factor(case_when(size == "<=20" ~ 0,
																														 size == "20-50" ~ 1,
																														 size == ">50" ~ 2)),
																	age = factor(age_q)) %>%
	select(dtime, death, size, age) %>%
	drop_na()

ggsurvplot(survfit(Surv(dtime, death) ~ size, data = rotterdam))
ggsurvplot(survfit(Surv(dtime, death) ~ age, data = rotterdam))

Dagum <- gamlss(Surv(dtime, death) ~ size + age, sigma.fo = ~ 1,
								data = rotterdam,
								family = cens(RD), method = RS(1000))
Gamma <- gamlss(Surv(dtime, death) ~ size + age, sigma.fo = ~ 1,
								data = rotterdam,
								family = cens(GA), method = RS(100))
Wei <- gamlss(Surv(dtime, death) ~ size + age, sigma.fo = ~ 1,
							data = rotterdam,
							family = cens(WEI3), method = RS(100))
Bs <- gamlss(Surv(dtime, death) ~ size + age, sigma.fo = ~ 1,
						 data = rotterdam,
						 family = cens(BS2), method = RS(100))
Ig <- gamlss(Surv(dtime, death) ~ size + age, sigma.fo = ~ 1,
						 data = rotterdam,
						 family = cens(IG), method = RS(100))

AICs <- list(Dagum = Dagum$aic, Gamma = Gamma$aic, Wei = Wei$aic, Bs = Bs$aic)
BICs <- list(Dagum = Dagum$sbc, Gamma = Gamma$sbc, Wei = Wei$sbc, Bs = Bs$sbc)
