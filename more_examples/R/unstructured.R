#Simple application of fitting mvn data with unstructured corr matrix

library(MASS);
library(RTMB);

#****** Data simulation
sd=diag(c(0.3,0.2,0.15));
cm=matrix(c(1.0 ,-0.5, -0.4,
               -0.5, 1.0, 0.4,
               -0.4, 0.4, 1.0), nrow=3,ncol=3)
Sigma=sd%*%cm%*%sd;

data=mvrnorm(n = 20, mu=c(0,0,0), 
        Sigma=Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE);
#***** End simulation code

# Problem set up;
us = unstructured(3); #to translate corr parameters to corr matrix
th=c(0,0,0);
mn=c(0,0,0);
logsd=c(log(.2),log(.2),log(.2));

dlst=list(xm=data);
plst=list(mn=mn,logsd=logsd,th=th);

f = function(par){
  getAll(dlst,par);
  s=exp(logsd);
  cm=us$corr(th);
  S=diag(s)%*%cm%*%diag(s);
  res = -sum(dmvnorm(x = xm, mu=mn,Sigma=S,log=TRUE));
  res
}
f(plst);

obj=MakeADFun(f,plst);
fit=nlminb(obj$par,obj$fn,obj$gr);
fit$convergence
sdr=sdreport(obj);
est=as.list(sdr,"Est");#get estimated pars as list
attr(est,"what")=NULL;
us$corr(est$th); #get and print the estimated corr matrix
