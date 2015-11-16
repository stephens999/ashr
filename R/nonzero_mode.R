# This file contains (experimental) methods for optimizing over a non-zero mode 

#' @title Estimate unimodal nonzero mean of a mixture model by EM algorithm
#'
#' @description Given the data, standard error of the data and standard deviations of the Gaussian mixture model, estimate the mean of a unimodal Gaussian mixture by an EM algorithm.
#'
#' @details Fits a k component mixture model \deqn{f(x|\pi) = \sum_k \pi_k f_k(x)} to independent
#' and identically distributed data \eqn{x_1,\dots,x_n}. 
#' Estimates unimodal mean \eqn{\mu} by EM algorithm. Uses the SQUAREM package to accelerate convergence of EM. Used by the ash main function; there is no need for a user to call this 
#' function separately, but it is exported for convenience. 
#' 
#' @param betahat a p vector of estimates.
#' @param sebetahat a p vector of corresponding standard errors.
#' @param mixsd vector of sds for underlying mixture components.
#' @param mixcompdist  distribution of components in mixture ( "uniform","halfuniform" or "normal").
#' @param df appropriate degrees of freedom for (t) distribution of betahat/sebetahat, default is NULL(Gaussian).
#' @param pi.init the initial value of \eqn{\pi} to use. If not specified defaults to (1/k,...,1/k).
#' @param control A list of control parameters for the SQUAREM algorithm, default value is set to be   control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE). 

#' 
#' @return A list, including the estimates (\eqn{\mu}) and (\eqn{\pi}), the log likelihood for each iteration (NQ)
#' and a flag to indicate convergence.
#' 
#' @export
#' 
#' 
nonzeromodeEM = function(betahat, sebetahat, mixsd, mixcompdist, df=NULL, pi.init=NULL,control=list()){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)  
  
  if(is.null(pi.init)){
    pi.init = rep(1/length(mixsd),length(mixsd))# Use as starting point for pi
  }
  else{
    pi.init=rgamma(length(mixsd),1,1)
  }
  
  if(controlinput$trace==TRUE){tic()}
  
  #tol in squarem needs special attention as it compares the difference of estimate fixed point (\mu),
  #thus, the absolute tol is not of our interest.here we set it to be relative tol, and our reference 
  #level would be the sample mean, the naive estimator,we keep significant figures to be 6 for fast 
  #convergence failing to set tol to adapt to data would result in keeping the fixed point iteration 
  #running till maxiter, and NaN would be produced, making the method inapplicable
  controlinput$tol=100*max(mean(betahat),1)*controlinput$tol
  
  if(length(sebetahat)==1){  sebetahat = rep(sebetahat,length(betahat))  }
  if(!is.element(mixcompdist,c("normal","uniform","halfuniform"))){
    stop("Error: invalid type of mixcompdist occcur in nonzeromodeEM()")
  }
  
  if(mixcompdist=="normal" & is.null(df)){
    g=normalmix(pi.init,rep(0,length(mixsd)),mixsd)
    mupi=c(mean(betahat),pi.init)
    res=squarem(par=mupi,fixptfn=nonzeromodeEMfixpoint,objfn=nonzeromodeEMobj,
                betahat=betahat,sebetahat=sebetahat,mixsd=mixsd,control=controlinput)
  }else if(mixcompdist=="normal" & !is.null(df)){
    stop("method comp_postsd of normal mixture not yet written for t likelihood")
  }else if(mixcompdist=="uniform"){
    g=unimix(pi.init,-mixsd,mixsd)
    mupi=c(mean(betahat),pi.init)    
    res=squarem(par=mupi,fixptfn=nonzeromodeEMoptimfixpoint,objfn=nonzeromodeEMoptimobj,
                betahat=betahat,sebetahat=sebetahat,g=g,df=df,control=controlinput)
  }else if(mixcompdist=="halfuniform"){
    g=unimix(c(pi.init, pi.init)/2,c(-mixsd,rep(0,length(mixsd))),c(rep(0,length(mixsd)),mixsd))
    mupi=c(mean(betahat),pi.init/2,pi.init/2)
    res=squarem(par=mupi,fixptfn=nonzeromodeEMoptimfixpoint,objfn=nonzeromodeEMoptimobj,
                betahat=betahat,sebetahat=sebetahat,g=g,df=df,control=controlinput)
  }
  
  if(controlinput$trace==TRUE){toc()}
  return(list(nonzeromode=res$par[1],pi=res$par[-1],NQ=-res$value.objfn,
              niter = res$iter, converged=res$convergence,post=res$par))
}


nonzeromodeEMoptimfixpoint = function(mupi,betahat,sebetahat,g,df){
  mu=mupi[1]
  pimean=normalize(pmax(0,mupi[-1])) #avoid occasional problems with negative pis due to rounding
  matrix_lik=t(compdens_conv(g,betahat-mu,sebetahat,df))
  m=t(pimean * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum=rowSums(m)
  classprob=m/m.rowsum #an n by k matrix
  pinew=normalize(colSums(classprob))
  munew=optimize(f=nonzeromodeEMoptim,interval=c(min(betahat),max(betahat)),pinew=pinew,betahat=betahat,sebetahat=sebetahat,g=g,df=df)$minimum
  mupi=c(munew,pinew)
  return(mupi)
}


nonzeromodeEMoptimobj = function(mupi,betahat,sebetahat,g,df){
  mu=mupi[1]
  pimean=normalize(pmax(0,mupi[-1])) #avoid occasional problems with negative pis due to rounding
  matrix_lik = t(compdens_conv(g,betahat-mu,sebetahat,df))
  m = t(pimean * t(matrix_lik))
  m.rowsum = rowSums(m)
  loglik = sum(log(m.rowsum))
  return(-loglik)
}


nonzeromodeEMoptim = function(mu,pinew,betahat,sebetahat,g,df){
  matrix_lik = t(compdens_conv(g,betahat-mu,sebetahat,df))
  pinew=normalize(pmax(0,pinew)) #avoid occasional problems with negative pis due to rounding
  m = t(pinew * t(matrix_lik))
  m.rowsum = rowSums(m)
  nloglik =- sum(log(m.rowsum))
  return(nloglik)
}


nonzeromodeEMfixpoint = function(mupi,betahat,sebetahat,mixsd){
  mu=mupi[1]
  pimean=normalize(pmax(0,mupi[-1])) #avoid occasional problems with negative pis due to rounding
  sdmat = sqrt(outer(sebetahat ^2,mixsd^2,"+")) 
  xmat=matrix(rep(betahat,length(mixsd)),ncol=length(mixsd))
  omegamatrix=t(t(dnorm(xmat,mean=mu,sd=sdmat))*pimean)
  omegamatrix=omegamatrix /rowSums(omegamatrix)
  pinew=normalize(colSums(omegamatrix))
  munew=sum(omegamatrix*xmat/(sdmat^2))/sum(omegamatrix/(sdmat^2))
  mupi=c(munew,pinew)
  return(mupi)
}

nonzeromodeEMobj = function(mupi,betahat,sebetahat,mixsd){
  mu=mupi[1]
  pimean=normalize(pmax(0,mupi[-1])) #avoid occasional problems with negative pis due to rounding
  sdmat = sqrt(outer(sebetahat ^2,mixsd^2,"+")) 
  xmat=matrix(rep(betahat,length(mixsd)),ncol=length(mixsd))
  omegamatrix=t(t(dnorm(xmat,mean=mu,sd=sdmat))*pimean)
  omegamatrix=omegamatrix /rowSums(omegamatrix)
  NegativeQ=-sum(omegamatrix*dnorm(xmat,mean=mu,sd=sdmat,log=TRUE))
  return(NegativeQ)
}

