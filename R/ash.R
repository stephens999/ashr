#TODO: Add nullcheck for VB?
#Separate out the optimization over sigma from the EM algorithm

#ash.repodir = scan(".ash.repodir.txt",what=character()) 
#source(file.path(ash.repodir, "/Rcode/mix.R"))
#source("mix.R")

#return the KL-divergence between 2 dirichlet distributions
#p,q are the vectors of dirichlet parameters of same lengths
diriKL = function(p,q){
  p.sum = sum(p)
  q.sum = sum(q)
  k = length(q)
  KL = lgamma(q.sum)-lgamma(p.sum)+sum((q-p)*(digamma(q)-digamma(rep(q.sum,k))))+sum(lgamma(p)-lgamma(q))
  return(KL)
}

#helper function for VBEM
VB.update = function(matrix_lik, pipost){
  avgpipost = matrix(exp(rep(digamma(pipost),n)-rep(digamma(sum(pipost)),k*n)),ncol=k,byrow=TRUE)
  classprob = avgpipost * matrix_lik
  classprob = classprob/rowSums(classprob) # n by k matrix
  B = sum(classprob*log(avgpipost*matrix_lik),na.rm=TRUE) - diriKL(prior,pipost) #negative free energy
  return(list(classprob=classprob,B=B))
}

#input matrix_lik is n by k matrix of p(obs n | comes from component k)
#prior: a k vector of dirichlet prior parameters
#output: post: a vector of the posterior dirichlet parameters
#B: the values of the likelihood lower bounds
#converged: boolean for whether it converged
#nullcheck: not implemented
VBEM = function(matrix_lik, prior, tol=0.0001, maxiter=5000){
  n=nrow(matrix_lik)
  k=ncol(matrix_lik)
  B = rep(0,maxiter)
  pipost = prior # Dirichlet posterior on pi
  
  avgpipost = matrix(exp(rep(digamma(pipost),n)-rep(digamma(sum(pipost)),k*n)),ncol=k,byrow=TRUE)
  classprob = avgpipost * matrix_lik
  classprob = classprob/rowSums(classprob) # n by k matrix  
  B[1] = sum(classprob*log(avgpipost*matrix_lik),na.rm=TRUE) - diriKL(prior,pipost) #negative free energy
 
  for(i in 2:maxiter){  
    pipost = colSums(classprob) + prior
    
    #Now re-estimate pipost
    avgpipost = matrix(exp(rep(digamma(pipost),n)-rep(digamma(sum(pipost)),k*n)),ncol=k,byrow=TRUE)
    classprob = avgpipost*matrix_lik
    classprob = classprob/rowSums(classprob) # n by k matrix
    
    B[i] = sum(classprob*log(avgpipost*matrix_lik),na.rm=TRUE) - diriKL(prior,pipost)
    
    if(abs(B[i]-B[i-1])<tol) break;
  }
  
  if(i>maxiter){i=maxiter}
   
  return(list(pihat = pipost/sum(pipost), B=B[1:i], niter = i, converged=(i<maxiter),post=pipost))
}
  

EM = function(matrix_lik, prior, pi.init = NULL,tol=0.0001, maxiter=5000,sigma.est=FALSE){
  n=nrow(matrix_lik)
  k=ncol(matrix_lik)
  B = rep(0,maxiter)
  pi = pi.init
  if(is.null(pi.init)){
    pi = rep(1/k,k)# Use as starting point for pi
  } 
  loglik = rep(0,maxiter)
  m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  loglik[1] = sum(log(m.rowsum))
  classprob = m/m.rowsum #an n by k matrix
  
  for(i in 2:maxiter){  
    pi = colSums(classprob) + prior-1
    pi = ifelse(pi<0,0,pi) #set any estimates that are less than zero, which can happen with prior<1, to 0
    pi = normalize(pi)
    
    #estimate sigma
    if(sigma.est==TRUE){
      sigmamin=min(sigmaavec)
      sigmamax=max(sigmaavec)
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
      matrix_lik = matrix_dens(betahat,sebetahat,sigmaavec)
    }
    
    #Now re-estimate pi
    m  = t(pi * t(matrix_lik)) 
    m.rowsum = rowSums(m)
    loglik[i] = sum(log(m.rowsum))
    classprob = m/m.rowsum
    
    if(abs(loglik[i]-loglik[i-1])<tol) break;
  }

  return(list(pihat = pi, B=loglik[1:i], niter = i, converged=(i<maxiter)))
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

EMest = function(betahat,sebetahat,g,sigma.est=FALSE,nullcheck=TRUE,prior=NULL,nc=NULL,VB=FALSE,ltol=0.0001, maxiter=5000){ 
 
  pi.init = g$pi
  k=ncomp(g)
  n = length(betahat)

  null.comp = 1; #which.min(sigmaavec) #which component is the "null"
  if(is.null(prior)){ # set up prior to be 1,1/(k-1),...,1/(k-1) to favour "null"
    prior = rep(1/(k-1),k)
    prior[null.comp] = 1
  }else if(prior=="uniform"){
    prior = rep(1,k)
  }
  
  matrix_lik = t(compdens_conv(g,betahat,sebetahat))
    
  if(VB==TRUE){
    EMfit=VBEM(matrix_lik,prior,ltol, maxiter)}
  else{
    EMfit = EM(matrix_lik,prior,pi.init,ltol, maxiter,sigma.est)
  }
  
  pi = EMfit$pihat     
  loglik = EMfit$B # actually return log lower bound not log-likelihood! 
  converged = EMfit$converged
  niter = EMfit$niter
  loglik.final = EMfit$B[niter]
  
  null.loglik = sum(log(matrix_lik[,null.comp]))  
  
  if(nullcheck==TRUE & VB==FALSE){ #null check doesn't work with VB yet
    if(null.loglik > loglik.final){ #check whether exceeded "null" likelihood where everything is null
      pi=rep(0,k)
      pi[null.comp]=1
      m  = t(pi * t(matrix_lik)) 
      m.rowsum = rowSums(m)
      loglik[niter] = sum(log(m.rowsum))
    }
  }
  
  g$pi=pi
  
  return(list(loglik=loglik[1:niter],null.loglik=null.loglik,
            matrix_lik=matrix_lik,converged = converged,g=g))
}



normalize = function(x){return(x/sum(x))}

#return the posterior on beta given a prior
#that is a mixture of normals (g)
#and observation betahat \sim N(beta,sebetahat)
#current matrix_lik is only for mu0=0, so would need to
#generalize that for general application
#INPUT: g: a normalmix with components indicating the prior
#       data, betahat (n vector), sebetahat (n vector)
#OUTPUT list (pi1,mu1,sigma1) whose components are each k by n matrices
#k is number of mixture components, n is number of observations
posterior_dist = function(g,betahat,sebetahat){
  pi0 = g$pi
  mu0 = g$mean
  sigma0 = g$sd
  
  k= length(pi0)
  n= length(betahat)
  
  pi1 = pi0 * t(matrix_dens(betahat,sebetahat,sigma0))
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


#return the "effective" estimate
#that is the effect size betanew whose z score betanew/se
#would give the same p value as betahat/se compared to a t with df
effective.effect=function(betahat,se,df){
  p = pt(betahat/se,df)
  qnorm(p,sd=se)
}

get_loglik = function(z.ash){
  return(rev(z.ash$fit$loglik)[1])
}

#compute corresponding q values from a vector of local fdr estimates
#INPUT: localfdr a vector of local fdr estimates
#OUTPUT: qvalue, a vector of q value estimates
qval.from.localfdr = function(localfdr){
  o = order(localfdr)
  qvalue=rep(NA,length(localfdr))
  qvalue[o] = (cumsum(sort(localfdr))/(1:sum(!is.na(localfdr))))
  return(qvalue)
}

# try to select a default range for the sigmaa values
# that should be used, based on the values of betahat and sebetahat
autoselect.sigmaavec = function(betahat,sebetahat){
  sigmaamin = min(sebetahat)/10 #so that the minimum is small compared with measurement precision
  if(all(betahat^2<sebetahat^2)){
    sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen, just arbitrarily use 4 normals.
  } else {
    sigmaamax = 2*sqrt(max(betahat^2-sebetahat^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2   
  }
  npoint = ceiling(log2(sigmaamax/sigmaamin))
  return(2^((-npoint):0) * sigmaamax)
}

#' @title Main Adaptive SHrinkage function
#'
#' @description Takes vectors of estimates (betahat) and their standard errors (sebetahat), and applies
#' shrinkage to them, using Empirical Bayes methods, to compute shrunk estimates for beta.
#'
#' @details See readme for more details
#' 
#' @param betahat, a p vector of estimates 
#' @param sebetahat, a p vector of corresponding standard errors
#' @param mixcompdist: distribution of components in mixture ("normal", "uniform" or "halfuniform")
#' @param df: appropriate degrees of freedom for (t) distribution of betahat/sebetahat
#' @param randomstart: bool, indicating whether to initialize EM randomly
#' @param usePointMass: bool, indicating whether to use a point mass at zero as one of components for a mixture distribution
#' @param onlylogLR (= FALSE) : bool, indicating whether to use this function to get logLR. Skip posterior prob, posterior mean, localfdr...
#' @param localfdr (=TRUE) : bool,  indicating whether to compute localfdr and q-value
#' @param auto (=TRUE): bool, whether to try to select the sigmaavec vector automatically (beta functionality)
#' @param sigma.est: bool, whether to estimate sigma rather than fixing (Beta version)
#' @param nc: number of components to use (only relevant when sigma.est=TRUE)
#' 
#' @return a list with elements fitted.g is fitted mixture
#' logLR : logP(D|mle(pi)) - logP(D|null)
#' 
#' @export
#' 
#' @examples 
#' beta = c(rep(0,100),rnorm(100))
#' sebetahat = abs(rnorm(200,0,1))
#' betahat = rnorm(200,beta,sebetahat)
#' beta.ash = ash(betahat, sebetahat)
#' summary(beta.ash)
#' plot(betahat,beta.ash$PosteriorMean,xlim=c(-4,4),ylim=c(-4,4))
#' 
#' 
#Things to do:
# check sampling routine
# check number of iterations
ash = function(betahat,sebetahat,mixcompdist = "normal",nullcheck=TRUE,df=NULL,randomstart=FALSE, usePointMass = FALSE, onlylogLR = FALSE, localfdr = TRUE, localfsr = TRUE, prior=NULL, sigmaavec=NULL, auto=TRUE, sigma.est=FALSE, nc=NULL, VB=FALSE){

  #if df specified, convert betahat so that bethata/sebetahat gives the same p value
  #from a z test as the original effects would give under a t test with df=df
  if(!is.null(df)){
    betahat = effective.effect(betahat,sebetahat,df)
  }	
    
  if(length(sebetahat)==1){
    sebetahat = rep(sebetahat,length(betahat))
  }
  if(length(sebetahat) != length(betahat)){
    stop("Error: sebetahat must have length 1, or same length as betahat")
  }
  if(is.null(sigmaavec)){
    sigmaavec = c(0.00025,0.0005,0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128,0.256,0.512,1.024,2.048,4.096,8.192) 
  }
  if(sigma.est==TRUE&!is.null(nc)){
    sigmaavec=2^(seq(-15,3,length.out=nc))
  }
  
  completeobs = (!is.na(betahat) & !is.na(sebetahat))
 if(sum(completeobs==0)){
    if(onlylogLR){
      return(list(pi=NULL, logLR = 0))
    }
    else{
      stop("Error: all input values are missing")
    }
  }  
if(auto==TRUE){
    sigmaavec= autoselect.sigmaavec(betahat[completeobs],sebetahat[completeobs])
  }
  if(usePointMass){
        sigmaavec = c(0,sigmaavec)
  }
  
  k=length(sigmaavec)  
  pi = rep(1, k)
  pi[1]=k
  pi=normalize(pi)
  if(randomstart){pi=rgamma(k,1,1)}
  
  if(!is.element(mixcompdist,c("normal","uniform","halfuniform"))) stop("Error: invalid type of mixcompdist")
  if(mixcompdist=="normal") g=normalmix(pi,rep(0,k),sigmaavec)
  if(mixcompdist=="uniform") g=unimix(pi,-sigmaavec,sigmaavec)
  if(mixcompdist=="halfuniform") g=unimix(c(pi,pi),c(-sigmaavec,rep(0,k)),c(rep(0,k),sigmaavec))
  
  pi.fit=EMest(betahat[completeobs],sebetahat[completeobs],g,sigma.est=sigma.est,prior=prior,nullcheck=nullcheck,nc=nc,VB=VB)  

  if(sigma.est==TRUE){
    sigmaavec=comp_sd(pi.fit$g)
  }
  
  if(onlylogLR){
	logLR = tail(pi.fit$loglik,1) - pi.fit$null.loglik
	return(list(pi=pi.fit$pi, logLR = logLR))
  }else{
#   	post = posterior_dist(pi.fit$g,betahat,sebetahat)
    n=length(betahat)
   	ZeroProb = rep(0,length=n)
    NegativeProb = rep(0,length=n)
    PosteriorMean = rep(0,length=n)
    PosteriorSD=rep(0,length=n)
    
   	ZeroProb[completeobs] = colSums(comppostprob(pi.fit$g,betahat[completeobs],sebetahat[completeobs])[comp_sd(pi.fit$g)==0,,drop=FALSE])     
   	NegativeProb[completeobs] = cdf_post(pi.fit$g, 0, betahat[completeobs],sebetahat[completeobs]) - ZeroProb[completeobs]
    PosteriorMean[completeobs] = postmean(pi.fit$g,betahat[completeobs],sebetahat[completeobs])
    PosteriorSD[completeobs] =postsd(pi.fit$g,betahat[completeobs],sebetahat[completeobs]) 
    
    #FOR MISSING OBSERVATIONS, USE THE PRIOR INSTEAD OF THE POSTERIOR
    ZeroProb[!completeobs] = sum(mixprop(pi.fit$g)[comp_sd(pi.fit$g)==0])
    NegativeProb[!completeobs] = mixcdf(pi.fit$g,0) 
    PosteriorMean[!completeobs] = mixmean(pi.fit$g)
    PosteriorSD[!completeobs] =mixsd(pi.fit$g)  
    PositiveProb =  1- NegativeProb-ZeroProb    
     
    
     
  	if(localfsr & localfdr){
   		localfsr = ifelse(PositiveProb<NegativeProb,PositiveProb+ZeroProb,NegativeProb+ZeroProb)
   		localfdr = 2* localfsr
      qvalue = qval.from.localfdr(localfdr)
  	}else{
   		localfdr=NULL
   		qvalue=NULL
  	}
   
    result = list(fitted.g=pi.fit$g,PosteriorMean = PosteriorMean,PosteriorSD=PosteriorSD,PositiveProb =PositiveProb,NegativeProb=NegativeProb, ZeroProb=ZeroProb,localfsr = localfsr, localfdr=localfdr,qvalue=qvalue,fit=pi.fit)
	  class(result)= "ash"
    return(result)

  }
  #if(nsamp>0){
  #  sample = posterior_sample(post,nsamp)
  #}
}

#' @title Summary method for ash object
#'
#' @description Print summary of fitted ash object
#'
#' @details See readme for more details
#' 
#' @export
#' 
summary.ash=function(a){
  print(a$fitted.g)
  print(tail(a$fit$loglik,1),digits=10)
  print(a$fit$converged)
}

#' @title Print method for ash object
#'
#' @description Print the fitted distribution of beta values in the EB hierarchical model
#'
#' @details None
#' 
#' @export
#' 
print.ash =function(a){
  print(a$fitted.g)
}

#' @title Plot method for ash object
#'
#' @description Plot the density of the underlying fitted distribution
#'
#' @details None
#' 
#' @export
#' 
plot.ash = function(a,xmin,xmax){
  x = seq(xmin,xmax,length=1000)
  y = density(a,x)
  plot(x,y,type="l")
}

#compute the predictive density of an observation
#given the fitted ash object a and the vector se of standard errors
#not implemented yet
predictive=function(a,se){
  
}


#' @title Plot method for ash object
#'
#' @description Return the density of the underlying fitted distribution
#'
#' @param a the fitted ash object
#' @param x the vector of locations at which density is to be computed
#'
#' @details None
#' 
#' @export
#' 
#'
density.ash=function(a,x){dens(a$fitted.g,x)}

#' @title Plot method for ash object
#'
#' @description Computed the cdf of the underlying fitted distribution
#'
#' @param a the fitted ash object
#' @param x the vector of locations at which cdf is to be computed
#' @param lower.tail (default=TRUE) whether to compute the lower or upper tail
#'
#' @details None
#' 
#' @export
#' 
#'
cdf.ash=function(a,x,lower.tail=TRUE){
 return(mixcdf(a$fitted.g,x,lower.tail))
}

