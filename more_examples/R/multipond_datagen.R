#This file includes code to generate and save simulated data as example
# to be used for a hierarchical model
#vonB parameterized with Linf, K, L2 rather than t0)
#"pond" specific vonB parameters logLinf, logK and LogL2 MVN

#simulated variation about pond specific vonB is normal
# with logCV = intercept + slope*(expected length at age)

library(MASS);

#Mean over ponds vonB params
logLinf_bar=log(500);
logK_bar=log(.35);
logL2_bar=log(175);

#var-cov matrix Sigma for vonB params (variation over ponds)
SD=diag(c(0.3,0.2,0.15));
corvb=matrix(c(1.0 ,-0.5, -0.4,
              -0.5, 1.0, 0.4,
              -0.4, 0.4, 1.0), nrow=3,ncol=3)
Sigma=SD%*%corvb%*%SD;

#generate pond-specific sample sizes, Z mortality rate, vonB params
nponds=20;
nfish=round(runif(n=nponds,min=100,max=500));
Z=rlnorm(n=nponds,meanlog=log(.45), sdlog=.2);
vonb_pond=mvrnorm(n = 20, mu=c(logLinf_bar,logK_bar,logL2_bar), 
                  Sigma=Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE);

#log CV for among fish variation in length is linear function of expect L given age
CV_int=log(0.4);
CV_slp=-(log(0.4)-log(0.05))/500;
#Note alt version just uses SDL
#SDL=20;


res=list();
for (p in 1:20){
 #age specific n Poisson with expectation based on pond-specific expected 
 # sample size (nfish) and pond-specific Z
 Na = rpois(n=length(2:25),nfish[p]*exp(-Z[p]*(2:25))/sum(exp(-Z[p]*(2:25))));
 #crude emulation of capping N by size category to avoid oversampling ages
 Na =round(pmin(Na,runif(length(Na),min=6.5,max=13.5)));
 Linf=exp(vonb_pond[p,1]);
 K=exp(vonb_pond[p,2]);
 L2=exp(vonb_pond[p,3]);
 tmplst=list();
 for (a in 2:25){
   mnL=L2+(Linf-L2)*(1-exp(-K*(a-2)));
   sd=exp(CV_int+CV_slp*mnL)*mnL; #log CV linear function of predicted length
  #  sd=SDL;
   if(Na[a-1]>0) {
      tmplst[[a-1]]=cbind(PID=rep(p),age=a,L=round(rnorm(n=Na[a-1],mean=mnL,sd=sd)));
   }
 }
   res=c(res,tmplst);
}   
     
 res2=data.frame(do.call(rbind,res));

 #Data visualization
 
 plot (res2$L~res2$age);
 
 for(i in levels(as.factor(res2$PID))) {
    plot(res2$L[res2$PID==i]~res2$age[res2$PID==i],main=paste0("PID= ",i));  
 }
 
 write.csv(x=res2,file="newmulti.csv",row.names=FALSE);
 
 
 
 
 
 
