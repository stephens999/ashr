#return ABF for vector of betahat and standard errors
ABF = function(betahat, sebetahat,sigmaa){
  T = betahat/sebetahat
  lambda = sebetahat^2/(sebetahat^2+sigmaa^2)
  return((sqrt(lambda) * exp(0.5*T^2 *(1-lambda))))
}

logABF = function(betahat,sebetahat,sigmaa){
  T = betahat/sebetahat
  lambda = sebetahat^2/(sebetahat^2+sigmaa^2)
  return(0.5*log(lambda) + 0.5*T^2 *(1-lambda))
  # the following line is same up to a constant, and probably faster:
  # return(dnorm(betahat,0,sqrt(sebetahat^2+sigmaa^2),log=TRUE))
}

#compute normal density for n-vector x
#at each of k values of mu and sigma
#OUTPUT: k by n matrix of normal densities
matdnorm = function (x, mu, sigma) 
{
  k=length(mu)
  n=length(x)
  d = matrix(rep(x,rep(k,n)),nrow=k)
  return(matrix(dnorm(d, mu, sigma),nrow=k))
}

#compute density of mixture of k normals
#for n-vector x
#INPUT: pi, mu and sigma are k vectors; x is n-vector
#OUTPUT: n-vector of densities
mixdnorm = function (x, pi, mu, sigma) 
{
  return (pi %*% matdnorm(x,mu,sigma))
}


#compute mixture log likelihood for data x, for mixture of k normals
#INPUT: x an n vector of data
#mixture component parameters pi, mu, sigma, each being k vectors
#
mixLoglik = function(x,pi,mu,sigma){
  dd=matdnorm(x,mu,sigma) 
  return(sum(log(pi %*% dd)))
}

#compute mixture log likelihood for data x, for mixture of k normals
#plus observation-specific errors with observation-specific variances
#ie x_s is sum_k pi_k N(x_s; mu_k, sigma^2_k + sigma^2_s)
# [optionally, use sigma^2_k * sigma^2_s if FUN="*"]
#INPUT: x an n vector of data
#mixture component parameters pi, mu, sigma, each being k vectors
#se, the observation-specific standard errrors
mixseLoglik = function(x,pi,mu,sigma,se,FUN="+"){
  k=length(mu)
  n=length(x)
  d = matrix(rep(x,rep(k,n)),nrow=k)
  s = sqrt(outer(sigma^2,se^2, FUN)) # k by n matrix of standard errors
  dd = matrix(dnorm(d, mu, s),nrow=k) 
  return(sum(log(pi %*% dd)))
}


#return matrix of ABFs for vector of sigma-a values
#normalized by maximum of each column
#betahat is n vector, sebetahat is n vector, sigmaavec is k vector
#return is n by k matrix of ABFs
matrixABF = function(betahat, sebetahat, sigmaavec){
  k = length(sigmaavec)
  n = length(betahat)
  labf = matrix(0,nrow=n, ncol=k)
  for(i in 1:k){
    labf[,i] = logABF(betahat,sebetahat,sigmaavec[i])
  }
  maxlabf = apply(labf, 1, max)
  labf = labf - maxlabf
  return(exp(labf))
}




#input abf is n by k matrix of p(obs n | comes from component k)
#prior: a k vector of dirichlet prior parameters
#output: post: a vector of the posterior dirichlet parameters
#B: the values of the likelihood lower bounds
#converged: boolean for whether it converged
#nullcheck: not implemented
oldVBEM = function(abf, prior, tol=0.0001, maxiter=5000,nullcheck=TRUE){
  n=nrow(abf)
  k=ncol(abf)
  B = rep(0,maxiter)
  pipost = prior # Dirichlet posterior on pi
  
  avgpipost = matrix(exp(rep(digamma(pipost),n)-rep(digamma(sum(pipost)),k*n)),ncol=k,byrow=TRUE)
  classprob = avgpipost * abf
  classprob = classprob/rowSums(classprob) # n by k matrix  
  B[1] = sum(classprob*log(avgpipost*abf)) - diriKL(prior,pipost) #negative free energy
  
  for(i in 2:maxiter){  
    pipost = colSums(classprob) + prior
    
    #Now re-estimate pipost
    avgpipost = matrix(exp(rep(digamma(pipost),n)-rep(digamma(sum(pipost)),k*n)),ncol=k,byrow=TRUE)
    classprob = avgpipost*abf
    classprob = classprob/rowSums(classprob) # n by k matrix
    
    B[i] = sum(classprob*log(avgpipost*abf)) - diriKL(prior,pipost)
    
    if(abs(B[i]-B[i-1])<tol) break;
  }
  
  if(i>maxiter){i=maxiter}
  
  return(list(post = pipost, B=B[1:i], niter = i, converged=(i<maxiter)))
}




#estimate mixture proportions of sigmaa by EM algorithm
#prior gives the parameter of a Dirichlet prior on pi
#(prior is used to encourage results towards smallest value of sigma when
#likelihood is flat)
#nullcheck indicates whether to check whether the loglike exceeds the null
#(may not want to use if prior is used)
#Introduced sigma.est with intention of implenting an option to esitmate
#sigma rather than fixing it, but this not yet implemented.
#VB provides an approach to estimate the approximate posterior distribution
#of mixture proportions of sigmaa by variational Bayes method
#(use Dirichlet prior and approximate Dirichlet posterior)

oldEMest = function(betahat,sebetahat,sigmaavec,pi,sigma.est=FALSE,nullcheck=TRUE,prior=NULL,nc=NULL,VB=FALSE,ltol=0.0001, maxiter=5000){ 
  if(!is.null(nc)&sigma.est==TRUE){
    sigmaavec=2^(seq(-15,3,length.out=nc))
  }
  k=length(sigmaavec)
  n = length(betahat)
  sigmamin=min(sigmaavec)
  sigmamax=max(sigmaavec)
  null.comp = which.min(sigmaavec) #which component is the "null"
  if(is.null(prior)){ # set up prior to be 1,1/(k-1),...,1/(k-1) to favour "null"
    prior = rep(1/(k-1),k)
    prior[null.comp] = 1
  }else if(prior=="uniform"){
    prior = rep(1,k)
  }
  
  abf = matrixABF(betahat,sebetahat,sigmaavec)
  
  if(VB==TRUE){
    vb=oldVBEM(abf,prior,ltol, maxiter)
    
    pi = vb$post/sum(vb$post) # use posterior mean to estimate pi    
    m  = t(pi * t(abf)) # abf is n by k; so this is also n by k
    m.rowsum = rowSums(m)
    classprob = m/m.rowsum
    loglik.final = sum(log(m.rowsum))
    null.loglik = sum(log(abf[,null.comp]))
    loglik = vb$B # actually return log lower bound not log-likelihood! 
    converged = vb$converged
    niter = vb$niter
  }else{
    loglik = rep(0,maxiter)
    m  = t(pi * t(abf)) # abf is n by k; so this is also n by k
    m.rowsum = rowSums(m)
    loglik[1] = sum(log(m.rowsum))
    classprob = m/m.rowsum #an n by k matrix
    
    for(i in 2:maxiter){  
      pi = colSums(classprob) + prior-1
      pi = ifelse(pi<0,0,pi) #set any estimates that are less than zero, which can happen with prior<1, to 0
      pi = normalize(pi)
      
      #estimate sigma
      if(sigma.est==TRUE){
        for(j in 1:k){
          pj=classprob[,j]
          f=function(x) sum((betahat^2/(sebetahat^2+x)^2-1/(sebetahat^2+x))*pj)
          if(f(sigmamin^2)<=0){
            sigmaavec[j]=sigmamin
          }else if(f(sigmamax^2)>=0){
            sigmaavec[j]=sigmamax
          }else{
            sigmaavec[j]=sqrt(uniroot(f,c(sigmamin^2,sigmamax^2))$root)          
          }
        }
        abf = matrixABF(betahat,sebetahat,sigmaavec)
      }
      
      #Now re-estimate pi
      m  = t(pi * t(abf)) 
      m.rowsum = rowSums(m)
      loglik[i] = sum(log(m.rowsum))
      classprob = m/m.rowsum
      
      if(abs(loglik[i]-loglik[i-1])<ltol) break;
    }
    null.loglik = sum(log(abf[,null.comp]))
    converged = (i< maxiter)
    niter= min(c(i,maxiter))
    loglik.final = loglik[niter]
    
    if(nullcheck==TRUE){
      if(null.loglik > loglik.final){ #check whether exceeded "null" likelihood where everything is null
        pi=rep(0,k)
        pi[null.comp]=1
        m  = t(pi * t(abf)) 
        m.rowsum = rowSums(m)
        loglik[niter] = sum(log(m.rowsum))
        classprob = m/m.rowsum
      }
    }
  }
  
  return(list(pi=pi,classprob=classprob,sigmaavec=sigmaavec,loglik=loglik[1:niter],null.loglik=null.loglik,
              abf=abf,converged = converged, temp1=sum(log(abf[,null.comp])),temp2=loglik.final))
}

normalize = function(x){return(x/sum(x))}

#return the posterior on beta given a prior
#that is a mixture of normals (pi0,mu0,sigma0)
#and observation betahat \sim N(beta,sebetahat)
#current ABF is only for mu0=0, so would need to
#generalize that for general application
#INPUT: priors: pi0, mu0, sigma0, all k vectors
#       data, betahat (n vector), sebetahat (n vector)
#OUTPUT list (pi1,mu1,sigma1) whose components are each k by n matrices
#k is number of mixture components, n is number of observations
old_posterior_dist = function(pi0,mu0,sigma0,betahat,sebetahat){
  k= length(pi0)
  n= length(betahat)
  
  pi1 = pi0 * t(matrixABF(betahat,sebetahat,sigma0))
  pi1 = apply(pi1, 2, normalize) #pi1 is now an k by n matrix
  
  #make k by n matrix versions of sigma0^2 and sebetahat^2
  # and mu0 and betahat
  s0m2 = matrix(sigma0^2,nrow=k,ncol=n,byrow=FALSE)
  sebm2 = matrix(sebetahat^2,nrow=k,ncol=n, byrow=TRUE)
  mu0m = matrix(mu0,nrow=k,ncol=n,byrow=FALSE)
  bhatm = matrix(betahat,nrow=k,ncol=n,byrow=TRUE)
  
  sigma1 = (s0m2*sebm2/(s0m2 + sebm2))^(0.5)  
  w = sebm2/(s0m2 + sebm2)
  mu1 = w*mu0m + (1-w)*bhatm
  
  #WHERE DATA ARE MISSING, SET POSTERIOR = PRIOR
  ismiss = (is.na(betahat) | is.na(sebetahat)) 
  pi1[,ismiss] = pi0
  mu1[,ismiss] = mu0
  sigma1[,ismiss] = sigma0
  
  return(list(pi=pi1,mu=mu1,sigma=sigma1))
}



#find mean and variance of a mixture of normals
#INPUT: x is a list with elements pi mu and sigma, each k by n matrices
#OUTPUT; the n vectors of mean and variances of mixtures correspondign to columns of pi, mu and sigma
normmix.mv=function(x){
  Ex = colSums(x$pi * x$mu)
  Ex2 = colSums(x$pi* (x$mu^2 + x$sigma^2))
  Vx = Ex2- Ex^2
  return(list(mean = Ex, var=Vx))
}  

#helper function for posterior_sample
#samples nsamp integers from 1:k according to a given prob vector
sample_component=function(p,nsamp){
  
  return(sample(length(p),nsamp,replace=TRUE,prob=p))
}

#m is a k by n matrix
#comp is a n vector of values in 1-k
#returns the comp[i]-th row of m[,i] 
extract_component=function(comp,m){
  return(m[cbind(comp,seq(comp))])
}

#returns matrix of nsamp samples from posterior
#computed using posterior_dist
# NOTE THIS IS UNTESTED, AND PROBABLY NOT WORKING YET...
posterior_sample = function(post,nsamp){
  component = as.vector(apply(post$pi,2,sample_component,nsamp=nsamp))
  k = ncol(post$pi)
  s = rep(1:k,rep(nsamp,k))
  index = cbind(component,s) #set up indices of mu and sigma to extract
  m = post$mu[index]
  ss = post$sigma[index]
  res = matrix(rnorm(length(m),mean=m,sd=ss),nrow=nsamp)
  return(res)
}


#find point estimates of beta from a posterior produced by posterior_dist
posterior_mean = function(post){
  return(colSums(post$pi * post$mu))
}

#return posterior of being <T (or >T) for a mixture of Gaussians
# each of pi1, mu1, sigma1 is a k by n matrix
# jth column provides parameters for jth mixture of gauusians 
# return an n vector of probabilities
pnormmix = function(T,pi1,mu1,sigma1,lower.tail=TRUE){
  return(apply(pi1 * pnorm(T,mu1,sigma1,lower.tail),2,sum))
}

#return posterior of being <T (or >T) for a mixture of Gaussians
# each of pi1, mu1, sigma1 is a k-vector
# jth element provides parameters for jth mixture of gauusians 
# return a probabilities
pnormmix.vec = function(T,pi1,mu1,sigma1,lower.tail=TRUE){
  return(sum(pi1 * pnorm(T,mu1,sigma1,lower.tail)))
}

#return the "effective" estimate
#that is the effect size betanew whose z score betanew/se
#would give the same p value as betahat/se compared to a t with df
effective.effect=function(betahat,se,df){
  p = pt(betahat/se,df)
  qnorm(p,sd=se)
}



#' @title Main Adaptive SHrinkage function (old version, for testing only)
#'
#' @description Takes vectors of estimates (betahat) and their standard errors (sebetahat), and applies
#' shrinkage to them, using Empirical Bayes methods, to compute shrunk estimates for beta.
#'
#' @details See readme for more details
#' 
#' @param betahat, a p vector of estimates 
#' @param sebetahat, a p vector of corresponding standard errors
#' @param mixcompdist: distribution of components in mixture ("normal", "uniform" or "halfuniform")
#' @param nullcheck: whether to check that any fitted model exceeds the "null" likelihood
#' in which all weight is on the first component

#' @param df: appropriate degrees of freedom for (t) distribution of betahat/sebetahat
#' @param randomstart: bool, indicating whether to initialize EM randomly
#' @param usePointMass: bool, indicating whether to use a point mass at zero as one of components for a mixture distribution
#' @param onlylogLR (= FALSE) : bool, indicating whether to use this function to get logLR. Skip posterior prob, posterior mean, localfdr...
#' @param localfdr (=TRUE) : bool,  indicating whether to compute localfdr, localfsr, and q-value
#' @param prior=NULL: vector of Dirichlet prior on mixture proportions (defaults to 1,1/k-1,...,1/k-1 to put more weight on first component)
#' @param sigmaavec=NULL: vector of parameters for underlying mixture components 
#' @param VB=FALSE: whether to use Variational Bayes to estimate mixture proportions (instead of EM to find MAP estimate)
#' 
#' 
#' @return a list with elements fitted.g is fitted mixture
#' logLR : logP(D|mle(pi)) - logP(D|null)
#' 
#' @export
#' 
#' 
#'
#main adaptive shrinkage function
#takes a vector of betahats and ses;
#fits a mixture of normals to it
# and returns posteriors
#INPUT: betahat (p vector); sebetahat (p vector of standard errors)
#df: degrees of freedome used to compute sebetahat
#randomstart: bool, indicating whether to initialize EM randomly
#usePointMass: bool, indicating whether to use a point mass at zero as one of components for a mixture distribution
#onlylogLR (= FALSE) : bool, indicating whether to use this function to get logLR. Skip posterior prob, posterior mean, localfdr...
#localfdr (=TRUE) : bool,  indicating whether to compute localfdr and q-value
#auto (=FALSE): bool, whether to try to select the sigmaavec vector automatically (beta functionality)
#OUTPUT: 
#logLR : logP(D|mle(pi)) - logP(D|null)
#Things to do: automate choice of sigmavec
# check sampling routin
# check number of iterations
oldash = function(betahat,sebetahat,nullcheck=TRUE,df=NULL,randomstart=FALSE, usePointMass = FALSE, onlylogLR = FALSE, localfdr = TRUE, localfsr = TRUE, prior=NULL, sigmaavec=NULL, auto=FALSE, sigma.est=FALSE, nc=NULL, VB=FALSE){
  
  #if df specified, convert betahat so that bethata/sebetahat gives the same p value
  #from a z test as the original effects would give under a t test with df=df
  if(!is.null(df)){
    betahat = effective.effect(betahat,sebetahat,df)
  }  
  
  if(length(sebetahat)==1){
    sebetahat = rep(sebetahat,length(betahat))
  }
  if(is.null(sigmaavec)){
    sigmaavec = c(0.00025,0.0005,0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128,0.256,0.512,1.024,2.048,4.096,8.192) 
  }
  
  completeobs = !is.na(betahat) & !is.na(sebetahat)
  if(auto==TRUE){
    sigmaavec= autoselect.mixsd(betahat[completeobs],sebetahat[completeobs])
  }
  if(usePointMass){
    sigmaavec = c(0,sigmaavec)
  }
  if(sigma.est==TRUE&!is.null(nc)){
    k=nc
  }else{
    k=length(sigmaavec)
  }
  pi = rep(1, k)
  pi[1]=k
  pi=normalize(pi)
  if(randomstart){pi=rgamma(k,1,1)}
  
  pi.fit=oldEMest(betahat[completeobs],sebetahat[completeobs],sigmaavec=sigmaavec,pi=pi,sigma.est=sigma.est,prior=prior,nullcheck=nullcheck,nc=nc,VB=VB)  
  if(sigma.est==TRUE){
    sigmaavec=pi.fit$sigmaavec
  }
  if(onlylogLR){
    logLR = pi.fit$temp2 - pi.fit$temp1
    return(list(pi=pi.fit$pi, logLR = logLR))
  }else{
    post = old_posterior_dist(pi.fit$pi,0,sigmaavec,betahat,sebetahat)
    PositiveProb = pnormmix(0,post$pi,post$mu,post$sigma,lower.tail=FALSE)
    ZeroProb = colSums(post$pi[sigmaavec==0,,drop=FALSE])
    NegativeProb =  1- PositiveProb-ZeroProb    
    PosteriorMean = posterior_mean(post)
    if(localfsr & localfdr){
      localfsr = ifelse(PositiveProb<NegativeProb,PositiveProb+ZeroProb,NegativeProb+ZeroProb)
      localfdr = 2* localfsr
      qvalue = qval.from.localfdr(localfdr)
    }else{
      localfdr=NULL
      qvalue=NULL
    }
    
    fitted.f= list(pi=pi.fit$pi,sigma=sigmaavec,mu=rep(0,k))
    result = list(post=post,fitted.f=fitted.f,PosteriorMean = PosteriorMean,PositiveProb =PositiveProb,NegativeProb=NegativeProb, ZeroProb=ZeroProb,localfsr = localfsr, localfdr=localfdr,qvalue=qvalue,fit=pi.fit)
    class(result)= "oldash"
    return(result)
    
  }
  #if(nsamp>0){
  #  sample = posterior_sample(post,nsamp)
  #}
}

print.oldash =function(a){
  a$fitted.f
}

plot.oldash = function(a,xmin,xmax){
  x = seq(xmin,xmax,length=1000)
  y = density(a,x)
  plot(x,y,type="l")
}

#compute the predictive density of an observation
#given the fitted ash object a and the vector se of standard errors
#not implemented yet
predictive=function(a,se){
  
}

#return the density of the fitted underlying hierarchical f
#a is an ash object
#x is a vector at which the density should be computed
density.oldash=function(a,x){mixdnorm(x,a$fitted.f$pi,a$fitted.f$mu,a$fitted.f$sigma)}

#return the cdf of the fitted underlying hierarchical f
#a is an ash object
#x is a vector at which the density should be computed
cdf.oldash=function(a,x,lower.tail=TRUE){
  return(vapply(x,pnormmix.vec, 0,pi1=a$fitted.f$pi,mu1=a$fitted.f$mu,sigma1=a$fitted.f$sigma,lower.tail=lower.tail))    
}

