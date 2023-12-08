# script for numerical search for the mean
# modified to include a grid search as well
# as search using nlminb without RTMB

#data
xvec=c(10.72,7.23,10.07,8.62,8.55);
#starting guesses for parameter
pars = c(8, log(sqrt(3)));

#My NLL function
NLL_fun = function(parms){
  -sum(dnorm(xvec,mean=parms[1],sd=exp(parms[2]),log=T))
}

#Produce vectors of 50 possible values of mu and sd
# over plausible ranges
mu_vals=seq(from=5,to=15, length.out=50)
sd_vals=seq(from=0.5, to=4, length.out=50);

# this is just left over code to make sure I can access
#  specific values of mu and sd and apply my NLL function
# before diving into a loop
i=5; j=7;
parms=c(mu_vals[i],log(sd_vals[j]));
tst=NLL_fun(parms);

#Create a 50 by 50 matrix and fill it with NLL values
# for each combo of mu and sd
NLL_Matrix=matrix(nrow=50,ncol=50);
for(i in 1:50){
  for(j in 1:50){
    parms=c(mu_vals[i],log(sd_vals[j]));
    NLL_Matrix[i,j]=NLL_fun(parms);
  }
}

#Finds the matrix row and col indices  
#  that had lowest negative log likelihood
index <- which(NLL_Matrix == min(NLL_Matrix), arr.ind = TRUE)
#these are the associated values of mu and sd that minimized NLL
mu_vals[index[1]];
sd_vals[index[2]];
#and the actual min NLL;
min(NLL_Matrix);

#Contour plot of NLL
contour(x=mu_vals,y=sd_vals,z=NLL_Matrix,levels=c(50,30,25,20,15,10,9))
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

