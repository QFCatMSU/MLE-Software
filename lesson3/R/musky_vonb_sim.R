# R script for fitting musky L@Age with RTMB
# profile CI and simulation add ons
# Adapted from admb and TMB classes by Bence and Brenden
set.seed(123456)
library(RTMB);

gmRdat = read.table("lesson3/data/musky_vonb.dat",head=T);

#Set up the data and starting value of parameters for RTMB
gmdat_real = list(len_obs=gmRdat[,"Length"],age=gmRdat[,"Age"]);
gmpar = list(log_linf=7,log_vbk=-1.6,t0=0,log_sd=4);

NLL_fun = function(par_lst){
  getAll(gmdat,par_lst);
  linf = exp(log_linf);
  vbk = exp(log_vbk);
  sd = exp(log_sd);
  len_obs=OBS(len_obs); # New line!
  len_pred = linf * (1 - exp(-vbk * (age - t0)));
  nll = -sum(dnorm(len_obs, len_pred, sd, TRUE));
  atage_pred = linf * (1 - exp(-vbk * ((1:11) - t0)))
  REPORT(atage_pred);
  nll
}

#Create object and fit model to original data
gmdat=gmdat_real;
obj=RTMB::MakeADFun(NLL_fun,gmpar);
fit=nlminb(obj$par, obj$fn, obj$gr);
fit$par;


simdat=gmdat_real;  #copy real data
#Simulate obs lengths and write to data copy
simdat$len_obs=obj$simulate()$len_obs;
gmdat=simdat;
objsim <- RTMB::MakeADFun(NLL_fun,gmpar);
simfit = nlminb(objsim$par, objsim$fn, objsim$gr);
simfit$par
fit$par
