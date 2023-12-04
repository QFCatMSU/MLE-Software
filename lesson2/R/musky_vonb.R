# R script for fitting musky L@Age with RTMB
# Adapted from admb and TMB classes by Bence and Brenden
library(RTMB);

gmRdat = read.table("lesson2/data/musky_vonb.dat",head=T);

#Set up the data and starting value of parameters for RTMB
gmdat = list(len_obs=gmRdat[,"Length"],age=gmRdat[,"Age"]);
gmpar = list(log_linf=7,log_vbk=-1.6,t0=0,log_sd=4);

#My NLL function
NLL_fun = function(par_lst){
  getAll(gmdat,par_lst);
  linf = exp(log_linf);
  vbk = exp(log_vbk);
  sd = exp(log_sd);
  len_pred = linf * (1 - exp(-vbk * (age - t0)));
  nll = -sum(dnorm(len_obs, len_pred, sd, TRUE));
  atage_pred = linf * (1 - exp(-vbk * ((1:11) - t0)))
  REPORT(atage_pred);
  nll
}

obj <- MakeADFun(NLL_fun,gmpar);

#Print report variables before fitting the model
GMreport=obj$report();
GMreport

## Call function minimizer
fit = nlminb(obj$par, obj$fn, obj$gr);

## Get parameter uncertainties and convergence diagnostics
sdr = sdreport(obj);
summary(sdr)

#Save the report variable list and print it after model fitting 
GMreport = obj$report();
GMreport;

#Do a likelihood profile CI for logLinf
#library(TMB);
#log_linf_prof = tmbprofile(obj,"log_linf");
#confint(log_linf_prof);
#plot(log_linf_prof)

