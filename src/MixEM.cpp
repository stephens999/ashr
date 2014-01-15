#include <Rcpp.h>
#include <cmath> //not sure this is necessary
#include <algorithm> 

using namespace Rcpp;

// [[Rcpp::export]]
List cxxMixEM(NumericMatrix matrix_lik, NumericVector prior, NumericVector pi_init, double tol=0.0001, int maxiter=5000){//note: no default pi_init=NULL  
  int n=matrix_lik.nrow(), k=matrix_lik.ncol(), j=0;
  bool converged=NA_LOGICAL;
  NumericVector pi(k);

  if(Rf_isNull(pi_init))
    std::fill(pi.begin(), pi.end(), 1./(double)k);
  else{
    pi=clone(pi_init);
    for (int i=0;i<k;i++)//set any estimates that are very small to be very small    
      pi[i]=std::max(1e-5, pi[i]); 
    pi=pi/sum(pi); //normalize pi
  }
  NumericMatrix m(n,k);
  NumericVector m_rowsum(n);
  NumericMatrix classprob(m);
  std::vector<double> loglik, priordens;
  loglik.reserve(maxiter);
  priordens.reserve(maxiter); 
  
  for (int i=0;i<k;i++){
    m.column(i)=pi[i]*matrix_lik.column(i);
    m_rowsum=m_rowsum+m.column(i);
  }
  for (int i=0;i<k;i++)//can vectorize this with sugar?
    classprob.column(i)=classprob.column(i)/m_rowsum;
  loglik.push_back(sum(log(m_rowsum)));
  priordens.push_back(sum((prior-1.)*log(pi)));

  for(j=1;j<maxiter;j++){
    //update pi
    //to do: can vectorize this with sugar?
    for (int i=0;i<k;i++)//set any estimates that are less than zero, which can happen with prior<1, to 0
      pi[i]=std::max(1e-5, sum(classprob.column(i))+prior[i]-1.);
    pi=pi/sum(pi); //normalize pi
    
    //Now re-estimate pi
    std::fill(m_rowsum.begin(), m_rowsum.end(), 0);
    for (int i=0;i<k;i++){
      m.column(i)=pi[i]*matrix_lik.column(i);
      m_rowsum=m_rowsum+m.column(i);
    }
    for (int i=0;i<k;i++)
      classprob.column(i)=classprob.column(i)/m_rowsum;
    loglik.push_back(sum(log(m_rowsum)));
    priordens.push_back(sum((prior-1.)*log(pi))); 
 
    converged=(bool) (std::abs(loglik[j]+priordens[j]-loglik[j-1]-priordens[j-1])<tol);
    if(converged)
      break;
  }
  
  return(List::create(Named("pihat")=pi,
                      Named("B")=loglik,
                      Named("niter")=wrap(min(j+1, maxiter)),
                      Named("converged")=wrap(converged)));
}
