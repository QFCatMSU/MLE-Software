# R script for fitting musky L@Age with RTMB
# Adapted from admb and TMB classes by Bence and Brenden
# Here we demo fitting by REML or ML and compare
# Show that you get same estimates except sd becomes approx unbiased estimator

library(RTMB)
gmdat <- read.table("lesson2/data/musky_vonb.dat", head = T)
# Set up the data and starting value of parameters for RTMB
dat <- list(len_obs = gmdat[, "Length"], age = gmdat[, "Age"])
par <- list(log_linf = 7, log_vbk = -1.6, t0 = 0, log_sd = 4)
# My NLL function
f <- function(par) {
  getAll(dat, par)
  len_obs <- OBS(len_obs)
  linf <- exp(log_linf)
  vbk <- exp(log_vbk)
  sd <- exp(log_sd)
  len_pred <- linf * (1 - exp(-vbk * (age - t0)))
  nll <- -sum(dnorm(len_obs, len_pred, sd, TRUE))
  atage_pred <- linf * (1 - exp(-vbk * ((1:11) - t0)))
  REPORT(atage_pred)
  nll
}

# Create two versions of "obj" one set up for ML estimation and one for REML
obj <- MakeADFun(f, par)
re <- c("log_linf", "log_vbk", "t0")
obj_reml <- MakeADFun(f, par, random = re)

# Call function minimizer
fit <- nlminb(obj$par, obj$fn, obj$gr)
fit_reml <- nlminb(obj_reml$par, obj_reml$fn, obj_reml$gr)

# Get parameter uncertainties and convergence diagnostics
sdr <- sdreport(obj)
summary(sdr)
sdr_reml <- sdreport(obj_reml)
summary(sdr_reml)

# For reml log_sd is now listed first in summary output
# - From RTMB perspective its the only fixed effect par!
# You won't see the other params if you just print sdr_reml
# sd is larger with reml (adjusting for negative bias of ML est)
# Notice that the other "parameters" have nearly identical values

# Numerical comparison of sds
sd <- exp(fit$par[4]) # ml sd
sd_reml <- exp(fit_reml$par[1]) # reml sd
n <- length(dat$len_obs)
# exclude sd (want pars that determine E(L$age))
p <- length(obj$par) - 1
# By known theory this should be same as REML
sd_adj <- sqrt((sd^2) * (n / (n - p)))
#print sd ests to compare
sd
sd_reml
sd_adj

