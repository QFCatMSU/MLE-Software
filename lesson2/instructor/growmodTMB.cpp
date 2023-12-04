#include <TMB.hpp>
#include <contrib/OSA_multivariate_dists-main/distr.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(A);
  DATA_VECTOR(L);
  DATA_VECTOR(agelst);
  
  PARAMETER(logLinf);
  PARAMETER(logK);
  PARAMETER(t0);
  PARAMETER(logSig);
  
  Type junk = 0;
  Type Linf = exp(logLinf);
  Type K = exp(logK);
  Type sig = exp(logSig);

  vector<Type> predL;
  vector<Type> resid;
  vector<Type> predlst;
  Type NLL;
    
  predL = Linf*(Type(1)-exp(-K*(A-t0)));
  resid = L-predL;
  
  predlst = Linf*(Type(1.)-exp(-K*(agelst-t0)));
  ADREPORT(predlst);

  NLL = -sum(dnorm(L,predL,sig,true));
  
  REPORT(Linf);
  REPORT(K);
  REPORT(sig);
  REPORT(predL);
  REPORT(resid);

  return NLL;
  
 
}
