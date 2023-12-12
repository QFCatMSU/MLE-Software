#Data generating  and estimating for age comp data

library(RTMB);

### code used to generate data
set.seed(123456);
log_rbar=1; #(arbitary)
Z=.15;
sd=0.1;
sel=c(.1,.15,.2,.3,.5,.68,.8,.85,.9,.94, .97, .99, .995, 1);
samp_size=c(100,100,10,50,50,50,10,200,200);
N=rep(0,14);
for(i in 1:14){
N[i]=exp(rnorm(n=1,mean=log_rbar-Z*(i-1),sd=sd));  
}
prop=N*sel/sum(N*sel);

samp_dat=matrix(nrow=14,ncol=9);
for (i in 1:9){
samp_dat[,i]=rmultinom(1,size=samp_size[i],prob=prop);
}

##### End data generation simulation code

real_dat_lst=list(samp_dat=samp_dat,samp_size=samp_size,sel=sel);
par_lst=list(Z=0.2);

NLL_fun=function(mypars){
  getAll(dat_lst,mypars);
  p=exp(-Z*(0:13))*sel; #Survival of recruit
  p=p/sum(p); #normalized to predict proportion at age
  NLL=0;
  for(i in 1:9){
    NLL=NLL-dmultinom(samp_dat[,i],prob=p,log=TRUE);
  }
  NLL;
}

dat_lst=real_dat_lst;
myobj = RTMB::MakeADFun(NLL_fun,par_lst);
opt <- nlminb(myobj$par,myobj$fn,myobj$gr);
sdrep = RTMB::sdreport(myobj);
summary(sdrep);

#Assume sel=1 for last four ages, monitonically increasing
logit=function(x){
  log(x/(1-x)); #x is on the real number line
}

inv_logit=function(x){ 
  1/(1+exp(-x)); # x is in (0,1)
}

adjust=rep(0.9,10); 
tr_adjust=logit(adjust);
par_lst=list(Z=Z,tr_adjust=tr_adjust);

NLL_fun=function(mypars){
  getAll(dat_lst,mypars);
  adj=inv_logit(tr_adjust);
  sel_use=rep(1,14);
  for(i in 10:1){sel_use[i]=sel_use[i+1]*adj[i]}
  ptmp=exp(-Z*(0:13))*sel_use; #Sel adj Survival of recruit
  p=ptmp/sum(ptmp); #normalized to predict proportion at age
  NLL=0;
  for(i in 1:9){
    NLL=NLL-dmultinom(samp_dat[,i],prob=p,log=TRUE);
    REPORT(sel_use)
  }
  NLL;
}
myobj = RTMB::MakeADFun(NLL_fun,par_lst);
opt <- nlminb(myobj$par,myobj$fn,myobj$gr);
opt$convergence
sdrep <- RTMB::sdreport(myobj);
sdrep;
myobj$report()
opt$convergence

#logistic selectivity
logistic=function(x,inf,slp){
1/(1+exp(-slp*(x-inf)))
}
logistic(1:14,5,0.5)

log_logist_inf=log(5);
log_logist_slp=log(0.5);
par_lst=list(Z=Z,log_logist_inf=log_logist_inf,log_logist_slp
             =log_logist_slp);

NLL_fun=function(mypars){
  getAll(dat_lst,mypars);
  logist_slp=exp(log_logist_slp);
  logist_inf=exp(log_logist_inf);
  sel_use=logistic(1:14,logist_inf,logist_slp);
  ptmp=exp(-Z*(0:13))*sel_use; #Sel adj Survival of recruit
  p=ptmp/sum(ptmp); #normalized to predict proportion at age
  NLL=0;
  for(i in 1:9){
    NLL=NLL-dmultinom(samp_dat[,i],prob=p,log=TRUE);
    REPORT(sel_use)
  }
  NLL;
}
myobj = RTMB::MakeADFun(NLL_fun,par_lst);
opt <- nlminb(myobj$par,myobj$fn,myobj$gr);
opt$convergence
sdrep <- RTMB::sdreport(myobj);
sdrep;
myobj$report()

# Now for fun, assume a linear form of Dirichlet-multinomial dsn
# Check out e.g., Fisch et al CJFAS 79:1745–1764 (2022)for
# log likelihood

# Dirichlet-multinomial is a compound distribution which assumes 
# probabilities come from a Dirichlet and given those probs, 
# n from multinomial.
  
ddirmult=function(x,alpha,log=FALSE){
  n=sum(x);
  salpha=sum(alpha);
  lden=lgamma(salpha)+lgamma(n+1)-lgamma(n+salpha) +
    sum(lgamma(x+alpha))-sum(lgamma(alpha))-sum(lgamma(x+1))
  if(log) ret=lden else ret=exp(lden);
  ret
}
# test my dirichlet multinomial density
 N_obs=dat_lst$samp_dat[,1];
 ddirmult(N_obs,N_obs*.8+0.01,log=T);
 
 par_lst=list(Z=Z,tr_scale=logit(0.8));
 
 NLL_fun=function(mypars){
   getAll(dat_lst,mypars);
   ptmp=exp(-Z*(0:13))*sel; #Sel adj Survival of recruit
   p=ptmp/sum(ptmp); #normalized to predict proportion at age
   scale=inv_logit(tr_scale);
   NLL=0;
   for(i in 1:9){
     samp_size=sum(samp_dat[,i])
     alpha=samp_size*scale*p;
     NLL=NLL-ddirmult(samp_dat[,i],alpha=alpha,log=TRUE);
   }
   NLL;
 }
 myobj = RTMB::MakeADFun(NLL_fun,par_lst);
 opt <- nlminb(myobj$par,myobj$fn,myobj$gr);
 opt$convergence
 sdrep <- RTMB::sdreport(myobj);
 sdrep;

 
 #Make the sample sizes way big for info content
 dat_lst=real_dat_lst;
 dat_lst$samp_dat=dat_lst$samp_dat*10;
 dat_lst$samp_size=dat_lst$samp_size*10;
 NLL_fun=function(mypars){
   getAll(dat_lst,mypars);
   ptmp=exp(-Z*(0:13))*sel; #Sel adj Survival of recruit
   p=ptmp/sum(ptmp); #normalized to predict proportion at age
   scale=inv_logit(tr_scale);
   NLL=0;
   for(i in 1:9){
     samp_size=sum(samp_dat[,i]);
     alpha=samp_size*scale*p;
     NLL=NLL-ddirmult(samp_dat[,i],alpha=alpha,log=TRUE);
   }
   NLL;
 }
 
 myobj = RTMB::MakeADFun(NLL_fun,par_lst);
 opt <- nlminb(myobj$par,myobj$fn,myobj$gr);
 opt$convergence
 sdrep <- RTMB::sdreport(myobj);
 sdrep;
 
 #In case someone wants to simulate
 # would have to do this manually, i.e., not using obj$simulate()
 rdirmult=function(n,size,alpha){
   #not vectorized!
   res_lst=list();
   k=length(alpha);
   for (i in 1:n){
     y=stats::rgamma(k,alpha,1);
     x=y/sum(y); #x~dirichlet
     res_lst[[i]]=as.vector(stats::rmultinom(1,size=size,prob=x));
   }
   do.call(rbind,res_lst);
 }
 rdirmult(1,sum(N_obs),alpha=N_obs*.5);
 
#"Manual simulation
NLL_fun=function(mypars){
  getAll(dat_lst,mypars);
  if (sim_code) sim_lst=list();
  ptmp=exp(-Z*(0:13))*sel; #Sel adj Survival of recruit
  p=ptmp/sum(ptmp); #normalized to predict proportion at age
  scale=inv_logit(tr_scale);
  NLL=0;
  for(i in 1:9){
    samp_size=sum(samp_dat[,i]);
    alpha=samp_size*scale*p;
    NLL=NLL-ddirmult(samp_dat[,i],alpha=alpha,log=TRUE);
    if(sim_code) sim_lst[[i]]=as.vector(rdirmult(1,samp_size,alpha));
  }
  if(sim_code) {
    sim_dat=do.call(cbind,sim_lst);
    ret=sim_dat;
  }
  else ret=NLL;
  ret;
}

#Showing this still works to estimate when sim_code is FALSE
sim_code=FALSE;
dat_lst=real_dat_lst;
dat_lst$samp_dat=dat_lst$samp_dat*10;
dat_lst$samp_size=dat_lst$samp_size*10;
myobj = RTMB::MakeADFun(NLL_fun,par_lst);
opt <- nlminb(myobj$par,myobj$fn,myobj$gr);
opt$convergence
sdrep <- RTMB::sdreport(myobj);
sdrep;

# Now using NLL_fun to simulate

#get ests as list
ests=as.list(sdrep,"Est");
attr(ests,"what")=NULL;
sim_code=TRUE;
sim_samp_dat=NLL_fun(ests);
sim_samp_dat;

#Fit model to simulated data
dat_lst=real_dat_lst;
dat_lst$samp_dat=sim_samp_dat;
dat_lst$samp_size=colSums(sim_samp_dat);
sim_code=FALSE;
simobj = RTMB::MakeADFun(NLL_fun,par_lst);
simopt <- nlminb(simobj$par,simobj$fn,simobj$gr);
simopt$convergence
sim_sdrep <- RTMB::sdreport(simobj);
sim_sdrep;
 