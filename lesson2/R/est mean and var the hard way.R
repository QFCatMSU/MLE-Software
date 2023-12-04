# script for numerical search for the mean

#data
xvec=c(10.72,7.23,10.07,8.62,8.55);
#starting guesses for parameter
pars = c(8, log(sqrt(3)));

#My NLL function
NLL_fun = function(parms){
  -sum(dnorm(xvec,mean=parms[1],sd=exp(parms[2]),log=T))
}

# Fit model using default finite diff derivatives
res = nlminb(pars,NLL_fun);
res;
res$par;
exp(res$par[2])^2; #back transform to get variance

#calculate the analytical param estimates
mean(xvec);
var(xvec);

#var gives the unbiased rather than MLE est of var
#convert to MLE (see bias example from lecture)
n=length(xvec);
var(xvec)*(n-1)/n;
