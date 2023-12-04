# script for RTMB search for the mean

#load RTMB
library(RTMB);

#data
xvec=c(10.72,7.23,10.07,8.62,8.55);
dat_lst=list(xvec=xvec)

#starting guesses for parameter
par_lst=list(mu=8,logsd=log(sqrt(3)));


#My NLL function
NLL_fun = function(par_lst){
  getAll(dat_lst,par_lst)
  -sum(dnorm(xvec,mean=mu,sd=exp(logsd),log=T))
}


obj <- MakeADFun(NLL_fun,par_lst);
opt <- nlminb(obj$par,obj$fn,obj$gr);
sdrep <- sdreport(obj)
summary(sdrep)

