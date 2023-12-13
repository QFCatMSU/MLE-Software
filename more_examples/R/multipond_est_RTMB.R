 #TMB program to fit hierarchical vonB for "realistic" data
 # Created Feb 2021 for MLE Software course by JRB
 # Model with params Linf, K, L2 (rather than t0)
 # {logLinf, logK, logL2} MV Normal
 # Modified July-Dec 2023 for RTMB

library(RTMB);

mp_df<-read.csv(file="more_examples/newmulti.csv");
head(mp_df);
 
 #Data visualization
 
 plot(mp_df$l~mp_df$age);
 
 for(i in levels(as.factor(mp_dat$pid))) {
    plot(mp_df$l[mp_df$pid==i]~mp_df$age[mp_df$pid==i], main = paste0("pid = ",i));  
 }

mp_dat=list(pid=mp_df$pid,obs_l=mp_df$l,age=mp_df$age);
 #find max length for each pond (use as starting Linf)
 max_l<-aggregate(l ~ pid, data = mp_df, max);
 
 mp_pars = list(log_linf_bar=log(500), log_linf=log(max_l$l), 
               log_k_bar=log(.3), log_k=rep(log(.3),20),
               log_l2_bar=log(175), log_l2=rep(log(175),20), 
               log_vb_sd=log(c(0.1,0.1,0.1)),
               theta=c(0,0,0),log_cv_int=log(0.1), log_cv_slp=0);

us = unstructured(3); 

f = function(pars){
 getAll(mp_dat, pars);
 linf_bar=exp(log_linf_bar);
 k_bar=exp(log_k_bar);
 l2_bar=exp(log_l2_bar);
 mn_vonb = c(log_linf_bar,log_k_bar,log_l2_bar);
 
 linf=exp(log_linf);
 k=exp(log_k);
 l2=exp(log_l2);
 vb_sd=exp(log_vb_sd);
 nobs=length(obs_l);
 nponds=length(log_l2);
 
 pred_l=sd=obs_l*0; #create pred_l and sd same length as obs_l and initialize with zeros

 for(i in 1:nobs){
   #predL = l2 +pred increase from l2 (vonB model starting at l2 at age=2)
   pred_l[i]=l2[pid[i]]+(linf[pid[i]]-l2[pid[i]])*(1-exp(-k[pid[i]]*(age[i]-2)));
   #log CV that determines SD for obs length linear function of predicted length
   sd[i]=exp(log_cv_int+log_cv_slp*pred_l[i])*pred_l[i];
 }

 nll = -sum(dnorm(x=obs_l,mean=pred_l,sd=sd,log=TRUE));
 
 #convert correlation params to correlation matrix
   #us = unstructured(3); 
 cm = us$corr(theta);
  #calculate var-cov matrix
 Sigma=diag(vb_sd) %*% cm %*% diag(vb_sd);
 vonb_pars=cbind(log_linf,log_k,log_l2);
 nprand = -sum(dmvnorm(x = vonb_pars, 
                       mu=mn_vonb,Sigma=Sigma,log=TRUE));
 nll+nprand;
 }
 
 re<-c("log_linf","log_k","log_l2");
 mpobj <- RTMB::MakeADFun(f, mp_pars,random=re);
 mpobj$report();

 ## Call function minimizer
 mpfit <- nlminb(mpobj$par, mpobj$fn, mpobj$gr);
 mpsdr <- sdreport(mpobj)
 summary(mpsdr);



 


 
 
