#todo
#
#' @title Main Adaptive SHrinkage function
#'
#' @description Takes vectors of estimates (betahat) and their standard errors (sebetahat), and applies
#' shrinkage to them, using Empirical Bayes methods, to compute shrunk estimates for beta.
#'
#' @details See readme for more details
#' 
#' @param betahat, a p vector of estimates 
#' @param sebetahat, a p vector of corresponding standard errors
#' @param method: specifies how ash is to be run. Can be "shrinkage" (if main aim is shrinkage) or "fdr" (if main aim is to assess fdr or fsr)
#' This is simply a convenient way to specify certain combinations of parameters: "shrinkage" sets pointmass=FALSE and prior="uniform";
#' "fdr" sets pointmass=TRUE and prior="nullbiased".
#' @param mixcompdist: distribution of components in mixture ("normal", "uniform" or "halfuniform")
#'
#' @param lambda1: multiplicative "inflation factor" for standard errors (like Genomic Control)
#' @param lambda2: additive "inflation factor" for standard errors (like Genomic Control)
#' @param nullcheck: whether to check that any fitted model exceeds the "null" likelihood
#' in which all weight is on the first component
#' @param df: appropriate degrees of freedom for (t) distribution of betahat/sebetahat
#' @param randomstart: bool, indicating whether to initialize EM randomly. If FALSE, then initializes to prior mean (for EM algorithm) or prior (for VBEM)
#' @param pointmass: bool, indicating whether to use a point mass at zero as one of components for a mixture distribution
#' @param onlylogLR: bool, indicating whether to use this function to get logLR. Skip posterior prob, posterior mean, lfdr...
#' @param prior: string, or numeric vector indicating Dirichlet prior on mixture proportions (defaults to "uniform", or 1,1...,1; also can be "nullbiased" 1,1/k-1,...,1/k-1 to put more weight on first component)
#' @param mixsd: vector of sds for underlying mixture components 
#' @param VB: whether to use Variational Bayes to estimate mixture proportions (instead of EM to find MAP estimate)
#' @param gridmult: the multiplier by which the default grid values for mixsd differ by one another. (Smaller values produce finer grids)
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
ash = function(betahat,sebetahat,method = c("shrink","fdr"), 
               mixcompdist = c("normal","uniform","halfuniform"),
               lambda1=1,lambda2=0,nullcheck=TRUE,df=NULL,randomstart=FALSE, 
               pointmass = FALSE, 
               onlylogLR = FALSE, 
               prior=c("uniform","nullbiased"), mixsd=NULL, VB=FALSE,gridmult=sqrt(2)){
  
    
  #If method is supplied, use it to set up defaults; provide warning if these default values
  #are also specified by user
  if(!missing(method)){
    method = match.arg(method) 
    if(method=="shrink"){
      if(missing(prior)){
        prior = "uniform"
      } else {
        cat("Warning: specification of prior overrides default for method shrink")
      }
      if(missing(pointmass)){
        pointmass=FALSE
      } else {
        cat("Warning: specification of pointmass overrides default for method shrink")
      }
    }
  
    if(method=="fdr"){
      if(missing(prior)){
        prior = "nullbiased"
      } else {
        cat("Warning: specification of prior overrides default for method fdr")
      }
      if(missing(pointmass)){
        pointmass=TRUE
      } else {
        cat("Warning: specification of pointmass overrides default for method fdr")
      }
    }  
  }
  
    mixcompdist = match.arg(mixcompdist)
    if(mixcompdist=="uniform" & pointmass==TRUE){
      stop("point mass not yet implemented for uniform or half-uniform")
    }
    if(mixcompdist=="halfuniform" & pointmass==TRUE){
      stop("point mass not yet implemented for uniform or half-uniform")
    }
    if(!is.numeric(prior)){
      prior = match.arg(prior)
    }
  
  
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
  
  completeobs = (!is.na(betahat) & !is.na(sebetahat))
  if(sum(completeobs)==0){
    if(onlylogLR){
      return(list(pi=NULL, logLR = 0))
    }
    else{
      stop("Error: all input values are missing")
    }
  }  
  
  if(is.null(mixsd)){
    mixsd= autoselect.mixsd(betahat[completeobs],sebetahat[completeobs],gridmult)
  }
  if(pointmass){
    mixsd = c(0,mixsd)
  }
  
  k=length(mixsd)  
  null.comp = which.min(mixsd) #which component is the "null"
  
  if(!is.numeric(prior)){
    if(prior=="nullbiased"){ # set up prior to favour "null"
      prior = rep(1,k)
      prior[null.comp] = 10 #prior 10-1 in favour of null
    }else if(prior=="uniform"){
      prior = rep(1,k)
    }
  }
  
  if(length(prior)!=k | !is.numeric(prior)){
    stop("invalid prior specification")
  }
  
  pi = prior #default is to initialize pi at prior (mean)
  if(randomstart){pi=rgamma(k,1,1)}
  
  if(!is.element(mixcompdist,c("normal","uniform","halfuniform"))) stop("Error: invalid type of mixcompdist")
  if(mixcompdist=="normal") g=normalmix(pi,rep(0,k),mixsd)
  if(mixcompdist=="uniform") g=unimix(pi,-mixsd,mixsd)
  if(mixcompdist=="halfuniform") g=unimix(c(pi,pi),c(-mixsd,rep(0,k)),c(rep(0,k),mixsd))
  
  pi.fit=EMest(betahat[completeobs],lambda1*sebetahat[completeobs]+lambda2,g,prior,null.comp=null.comp,nullcheck=nullcheck,VB=VB)  
  
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
    
    lfsr = ifelse(PositiveProb<NegativeProb,PositiveProb+ZeroProb,NegativeProb+ZeroProb)
    lfsra = ifelse(PositiveProb<NegativeProb,2*PositiveProb+ZeroProb,2*NegativeProb+ZeroProb)    
    lfdr = ZeroProb
    qvalue = qval.from.lfdr(lfdr)
        
    result = list(fitted.g=pi.fit$g,logLR =tail(pi.fit$loglik,1) - pi.fit$null.loglik,PosteriorMean = PosteriorMean,PosteriorSD=PosteriorSD,PositiveProb =PositiveProb,NegativeProb=NegativeProb, ZeroProb=ZeroProb,lfsr = lfsr,lfsra=lfsra, lfdr=lfdr,qvalue=qvalue,fit=pi.fit,lambda1=lambda1,lambda2=lambda2,call=match.call(),data=list(betahat = betahat, sebetahat=sebetahat))
    class(result)= "ash"
    return(result)
    
  }
  #if(nsamp>0){
  #  sample = posterior_sample(post,nsamp)
  #}
}


#' @title Estimate posterior distribution on mixture proportions of a mixture model by a Variational Bayes EM algorithm
#'
#' @description Given the individual component likelihoods for a mixture model, estimates the posterior on 
#' the mixture proportions by an VBEM algorithm. Used by the ash main function; there is no need for a user to call this 
#' function separately, but it is exported for convenience.
#'
#' @details Fits a k component mixture model \deqn{f(x|\pi) = \sum_k \pi_k f_k(x)} to independent
#' and identically distributed data \eqn{x_1,\dots,x_n}. 
#' Estimates posterior on mixture proportions \eqn{\pi} by Variational Bayes, 
#' with a Dirichlet prior on \eqn{\pi}. 
#' Algorithm adapted from Bishop (2009), Pattern Recognition and Machine Learning, Chapter 10.
#' 
#' @param matrix_lik: a n by k matrix with (j,k)th element equal to \eqn{f_k(x_j)}.
#' @param prior: a k vector of the parameters of the Dirichlet prior on \eqn{\pi}. Recommended to be rep(1,k)
#' @param post.init: the initial value of the posterior parameters. If not specified defaults to the prior parameters.
#' @param tol: the tolerance for convergence of log-likelihood bound.
#' @param maxiter: the maximum number of iterations performed
#' 
#' @return A list, whose components include point estimates (pihat), 
#' the parameters of the fitted posterior on \eqn{\pi} (pipost),
#' the bound on the log likelihood for each iteration (B)
#' and a flag to indicate convergence (converged).
#'  
#' @export
#' 
mixVBEM = function(matrix_lik, prior, post.init=NULL, tol=0.0001, maxiter=5000){
  n=nrow(matrix_lik)
  k=ncol(matrix_lik)
  B = rep(0,maxiter)
  pipost = post.init
  if(is.null(post.init)){
    pipost = prior # Dirichlet posterior on pi
  }
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
   
  return(list(pihat = pipost/sum(pipost), B=B[1:i], niter = i, converged=(abs(B[i]-B[i-1])<tol),post=pipost))
}
  

#' @title Estimate mixture proportions of a mixture model by EM algorithm
#'
#' @description Given the individual component likelihoods for a mixture model, estimates the mixture proportions by an EM algorithm.
#'
#' @details Fits a k component mixture model \deqn{f(x|\pi) = \sum_k \pi_k f_k(x)} to independent
#' and identically distributed data \eqn{x_1,\dots,x_n}. 
#' Estimates mixture proportions \eqn{\pi} by maximum likelihood, or by maximum a posteriori (MAP) estimation for a Dirichlet prior on $\pi$ 
#' (if a prior is specified).  Used by the ash main function; there is no need for a user to call this 
#' function separately, but it is exported for convenience.
#'
#' 
#' @param matrix_lik, a n by k matrix with (j,k)th element equal to \eqn{f_k(x_j)}.
#' @param prior, a k vector of the parameters of the Dirichlet prior on \eqn{\pi}. Recommended to be rep(1,k)
#' @param pi.init, the initial value of \eqn{\pi} to use. If not specified defaults to (1/k,...,1/k).
#' @param tol, the tolerance for convergence of log-likelihood.
#' @param maxiter the maximum number of iterations performed
#' 
#' @return A list, including the estimates (pihat), the log likelihood for each interation (B)
#' and a flag to indicate convergence
#'  
#' @export
#' 
mixEM = function(matrix_lik, prior, pi.init = NULL,tol=0.0001, maxiter=5000){
  n=nrow(matrix_lik)
  k=ncol(matrix_lik)
  B = rep(0,maxiter)
  pi = pi.init
  if(is.null(pi.init)){
    pi = rep(1/k,k)# Use as starting point for pi
  } 
  pi = ifelse(pi<1e-5,1e-5,pi) #set any estimates that are very small to be very small
  pi = normalize(pi)
  
  loglik = rep(0,maxiter)
  priordens= rep(0,maxiter)
  m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  loglik[1] = sum(log(m.rowsum))
  priordens[1] = sum((prior-1)*log(pi)) 
  classprob = m/m.rowsum #an n by k matrix
  
  for(i in 2:maxiter){  
    pi = colSums(classprob) + prior-1
    pi = ifelse(pi<1e-5,1e-5,pi) #set any estimates that are less than zero, which can happen with prior<1, to 0
    pi = normalize(pi)
        
    #Now re-estimate pi
    m  = t(pi * t(matrix_lik)) 
    m.rowsum = rowSums(m)
    loglik[i] = sum(log(m.rowsum))
    priordens[i] = sum((prior-1)*log(pi)) 
    classprob = m/m.rowsum
    
    
    if(abs(loglik[i]+priordens[i]-loglik[i-1]-priordens[i-1])<tol) break;
  }

  return(list(pihat = pi, B=loglik[1:i], 
              niter = i, converged=(abs(loglik[i]+priordens[i]-loglik[i-1]-priordens[i-1])<tol)))
}



#estimate mixture proportions of sigmaa by EM algorithm
#prior gives the parameter of a Dirichlet prior on pi
#(prior is used to encourage results towards smallest value of sigma when
#likelihood is flat)
#nullcheck indicates whether to check whether the loglike exceeds the null
#(may not want to use if prior is used)
#VB provides an approach to estimate the approximate posterior distribution
#of mixture proportions of sigmaa by variational Bayes method
#(use Dirichlet prior and approximate Dirichlet posterior)

EMest = function(betahat,sebetahat,g,prior,null.comp=1,nullcheck=TRUE,VB=FALSE,ltol=0.0001, maxiter=5000){ 
 
  pi.init = g$pi
  k=ncomp(g)
  n = length(betahat)
  
  matrix_lik = t(compdens_conv(g,betahat,sebetahat))
    
  if(VB==TRUE){
    EMfit=mixVBEM(matrix_lik,prior,pi.init,ltol, maxiter)}
  else{
    EMfit = mixEM(matrix_lik,prior,pi.init,ltol, maxiter)
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


#' @title Estimate mixture proportions of a mixture model by EM algorithm
#'
#' @description Return the posterior on beta given a prior (g) that is a mixture of normals (class normalmix) 
#' and observation betahat \sim N(beta,sebetahat)
#'
#' @details This can be used to obt
#'
#' @param g: a normalmix with components indicating the prior; works only if g has means 0
#' @param betahat (n vector of observations) 
#' @param sebetahat (n vector of standard errors/deviations of observations)
#' 
#' @return A list, (pi1,mu1,sigma1) whose components are each k by n matrices
#' where k is number of mixture components in g, n is number of observations in betahat
#' 
#' @export
#' 
#' 
posterior_dist = function(g,betahat,sebetahat){
  if(class(g)!="normalmix"){
    stop("Error: posterior_dist implemented only for g of class normalmix")
  }
  pi0 = g$pi
  mu0 = g$mean
  sigma0 = g$sd  
  k= length(pi0)
  n= length(betahat)
  
  if(!all.equal(g$mean,rep(0,k))) stop("Error: posterior_dist currently only implemented for zero-centered priors")
  
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

#return matrix of densities of observations (betahat) 
# assuming betahat_j \sim N(0, sebetahat_j^2 + sigmaavec_k^2)
#normalized by maximum of each column
#INPUT
#betahat is n vector, 
#sebetahat is n vector, 
#sigmaavec is k vector
#return is n by k matrix of the normal likelihoods, 
# with (j,k)th element the density of N(betahat_j; mean=0, var = sebetahat_j^2 + sigmaavec_k^2)
#normalized to have maximum 1 in each column
matrix_dens = function(betahat, sebetahat, sigmaavec){
  k = length(sigmaavec)
  n = length(betahat)
  ldens = dnorm(betahat,0,sqrt(outer(sebetahat^2,sigmaavec^2,FUN="+")),log=TRUE)
  maxldens = apply(ldens, 1, max)
  ldens = ldens - maxldens
  return(exp(ldens))
}

#return the "effective" estimate
#that is the effect size betanew whose z score betanew/se
#would give the same p value as betahat/se compared to a t with df
effective.effect=function(betahat,se,df){
  p = pt(betahat/se,df)
  qnorm(p,sd=se)
}


#' @title Function to compute q values from local false discovery rates
#'
#' @description Computes q values from a vector of local fdr estimates
#'
#' @details The q value for a given lfdr is an estimate of the (tail) False Discovery Rate 
#' for all findings with a smaller lfdr, and is found by the average of the lfdr for
#' all more significant findings. See Storey (2003), Annals of Statistics, for definition of q value.  
#' 
#' 
#' @param lfdr, a vector of local fdr estimates
#'
#' @return vector of q values
#' 
#' @export
qval.from.lfdr = function(lfdr){
  o = order(lfdr)
  qvalue=rep(NA,length(lfdr))
  qvalue[o] = (cumsum(sort(lfdr))/(1:sum(!is.na(lfdr))))
  return(qvalue)
}

# try to select a default range for the sigmaa values
# that should be used, based on the values of betahat and sebetahat
# mult is the multiplier by which the sds differ across the grid
autoselect.mixsd = function(betahat,sebetahat,mult){
  sigmaamin = min(sebetahat)/10 #so that the minimum is small compared with measurement precision
  if(all(betahat^2<sebetahat^2)){
    sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
  } else {
    sigmaamax = 2*sqrt(max(betahat^2-sebetahat^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2   
  }
  npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
  return(mult^((-npoint):0) * sigmaamax)
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
plot.ash = function(a,xmin,xmax,...){
  x = seq(xmin,xmax,length=1000)
  y = density(a,x)
  plot(x,y,type="l",...)
}

#compute the predictive density of an observation
#given the fitted ash object a and the vector se of standard errors
#not implemented yet
predictive=function(a,se){
  
}


#' @title Get fitted loglikelihood for ash object
#'
#' @description Return the log-likelihood of the data under the fitted distribution
#'
#' @param a the fitted ash object
#'
#' @details None
#' 
#' @export
#' 
#'
get_loglik = function(a){
  return(tail(a$fit$loglik,1))
}

#' @title Get pi0 estimate for ash object
#'
#' @description Return estimate of the null proportion, pi0
#'
#' @param a the fitted ash object
#'
#' @details Extracts the estimate of the null proportion, pi0, from the object a
#' 
#' @export
#' 
get_pi0 = function(a){
  return(ifelse(comp_sd(a$fitted.g)[1]==0,a$fitted.g$pi[1],0))
}

#' @title Compute loglikelihood for data from ash fit
#'
#' @description Return the log-likelihood of the data betahat, with standard errors betahatsd, 
#' under the fitted distribution in the ash object
#'
#' @param a the fitted ash object
#' @param betahat the data
#' @param betahatsd the observed standard errors
#' 
#' @details None
#' 
#' @export
#' 
#'
loglik.ash = function(a,betahat,betahatsd){
  return(loglik_conv(a$fitted.g,betahat, betahatsd))
}

#' @title Density method for ash object
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
density.ash=function(a,x){list(x=x,y=dens(a$fitted.g,x))}

#' @title cdf method for ash object
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
 return(list(x=x,y=mixcdf(a$fitted.g,x,lower.tail)))
}


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
