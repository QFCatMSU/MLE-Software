# R script for fitting growmod example with TMB
# Dev Dec 2018, Jan 2019 by JRB for MLE Software class

library(TMB);
gmRdat = read.table("growmodTMB.dat",header=T);

#Set up the data and starting value of parameters for TMB
gmdat = list(L=gmRdat[,"Length"],A=gmRdat[,"Age"]);
gmpar = list(logLinf=7,logK=-1.6,t0=0,logSig=4);

#Revise data list for including agelst (don't use until you change the template for this)
agelst = 1:11;
gmdat = list(L=gmRdat[,"Length"],A=gmRdat[,"Age"],agelst=agelst);

#Map with logLinf not estimated -- do not run or reset logLinf to factor(1) if you want to estimate it!
gmmap = list(logLinf=factor(NA));

## Compile and load the model
compile("growmodTMB.cpp");
dyn.load(dynlib("growmodTMB"));

## Make a function object
objGM = MakeADFun(gmdat, gmpar, DLL="growmodTMB");

## Make a function object that uses a map.  Do not use unless you intend the map!
objGM = MakeADFun(gmdat, gmpar, map=gmmap, DLL="growmodTMB");

#Print report variables before fitting the model
objGM$report();

#Save the report variable list and print one before fitting the model
GMreport = objGM$report();
GMreport$Linf;

## Call function minimizer
fit = nlminb(objGM$par, objGM$fn, objGM$gr);

## Get parameter uncertainties and convergence diagnostics
sdr = sdreport(objGM);
sdr;

#Save the report variable list and print it after model fitting 
GMreport = objGM$report();
GMreport$Linf;

#Particularly useful when you have used ADREPORT
summary(sdr);

#Do a likelihood profile CI for logLinf
logLinfprof = tmbprofile(objGM,"logLinf");
confint(logLinfprof);

# Code for running with bounds
# create a function object with starting values
objGMb = MakeADFun(gmdat, gmpar, DLL="growmodTMB");
#create upper and lower bounds
lower = c(4,-5,-5,0);
upper = c(10,1,5,10);

#fit with bounds -- note use of objGMb
bfit = nlminb(objGMb$par, objGMb$fn, objGMb$gr,upper=upper, lower=lower);
sdrb = sdreport(objGMb);
summary(sdrb);