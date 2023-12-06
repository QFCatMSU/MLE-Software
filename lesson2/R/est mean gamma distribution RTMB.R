# script for est gamma distribution params by RTMB

#load RTMB
library(RTMB);

#data
xvec=c(10.72,7.23,10.07,8.62,8.55);
dat_lst=list(xvec=xvec)

#X~gamma(shape,scale) then (from help(gamma))
#E(X)=shape*scale V(X)=shape*scale^2
#scale=V(X)/E(X)
#shape=E(X)/scale
#starting guesses for parameter
par_lst=list(logscale=log(3/8),logshape= log(8/(3/8)));


#My NLL function
NLL_fun = function(par_lst){
  getAll(dat_lst,par_lst);
  shape=exp(logshape);
  scale=exp(logscale);
  logSD=log(sqrt(shape*scale^2));
  mu=shape*scale;
  ADREPORT(logSD);
  ADREPORT(mu);
  -sum(dgamma(xvec,scale=scale,shape=shape,log=T))
}


obj <- MakeADFun(NLL_fun,par_lst);
opt <- nlminb(obj$par,obj$fn,obj$gr);
sdrep <- sdreport(obj);
sumsd=summary(sdrep);
sdrep;
sumsd;
obj$env$value.best #NLL
#MLEs (same as in summary(sdrep))

#Pull out parameter estimates and standard errors
est=as.list(sdrep,"Est");
est;
se=as.list(sdrep,"Std");
se;
# strip out "what" attribute from est
attr(est,"what")=NULL;
est; #same form as the starting value list
