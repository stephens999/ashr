#' @useDynLib ashr
#todo
#' @import truncnorm SQUAREM doParallel pscl Rcpp
#
#
#' @title Main Adaptive Shrinkage function
#'
#' @description Takes vectors of estimates (betahat) and their standard errors (sebetahat), and applies
#' shrinkage to them, using Empirical Bayes methods, to compute shrunk estimates for beta.
#'
#' @details See readme for more details
#' 
#' @param betahat  a p vector of estimates 
#' @param sebetahat a p vector of corresponding standard errors
#' @param method specifies how ash is to be run. Can be "shrinkage" (if main aim is shrinkage) or "fdr" (if main aim is to assess fdr or fsr)
#' This is simply a convenient way to specify certain combinations of parameters: "shrinkage" sets pointmass=FALSE and prior="uniform";
#' "fdr" sets pointmass=TRUE and prior="nullbiased".
#' @param mixcompdist distribution of components in mixture ( "uniform","halfuniform" or "normal"), the default value would be "uniform"
#'
#' @param lambda1  multiplicative "inflation factor" for standard errors (like Genomic Control)
#' @param lambda2  additive "inflation factor" for standard errors (like Genomic Control)
#' @param nullcheck  whether to check that any fitted model exceeds the "null" likelihood
#' in which all weight is on the first component
#' @param df appropriate degrees of freedom for (t) distribution of betahat/sebetahat, default is NULL(Gaussian)
#' @param nullweight scalar, the weight put on the prior under "nullbiased" specification, see \code{prior}
#' @param randomstart logical, indicating whether to initialize EM randomly. If FALSE, then initializes to prior mean (for EM algorithm) or prior (for VBEM)
#' @param nonzeromode logical, indicating whether to use a non-zero unimodal mixture(default is "FALSE")
#' @param pointmass logical, indicating whether to use a point mass at zero as one of components for a mixture distribution
#' @param onlylogLR logical, indicating whether to use this function to get logLR. Skip posterior prob, posterior mean, lfdr...
#' @param prior string, or numeric vector indicating Dirichlet prior on mixture proportions (defaults to "uniform", or (1,1...,1); also can be "nullbiased" (nullweight,1,...,1) to put more weight on first component)
#' @param mixsd vector of sds for underlying mixture components 
#' @param VB whether to use Variational Bayes to estimate mixture proportions (instead of EM to find MAP estimate), see \code{\link{mixVBEM}} and \code{\link{mixEM}}
#' @param gridmult the multiplier by which the default grid values for mixsd differ by one another. (Smaller values produce finer grids)
#' @param minimaloutput if TRUE, just outputs the fitted g and the lfsr (useful for very big data sets where memory is an issue) 
#' @param multiseqoutput if TRUE, just outputs the fitted g, logLR, PosteriorMean, PosteriorSD, function call and df
#' @param g the prior distribution for beta (usually estimated from the data; this is used primarily in simulated data to do computations with the "true" g)
#' @param K An integer denoting the order of the SQUAREM scheme. Default is 1,i.e. first-order schemes, which is adequate for most problems. K=2,3 may provide greater speed in some problems, although they are less reliable than the first-order schemes.
#' @param maxiter maximum number of iterations of the EM algorithm.
#' @param retol the relelative precision for the mode of mixture when nonzeromode=TRUE, the default value is 1e-5.
#' @param trace a logical variable denoting whether some of the intermediate results of iterations should be displayed to the user. Default is FALSE.
#' @param cxx flag to indicate whether to use the c++ (Rcpp) version. After application of Squared extrapolation methods for accelerating fixed-point iterations (R Package "SQUAREM"), the c++ version is no longer faster than non-c++ version, thus we do not recommend using this one, and might be removed at any point. 
#' @param model c("EE","ES") specifies whether to assume exchangeable effects (EE) or exchangeable standardized effects (ES).
#' @param control A list of control parameters specifing any changes to default values of algorithm control parameters. Full names of control list elements must be specified, otherwise, user-specifications are ignored. See Details.
#' 
#' @return ash returns an object of \code{\link[base]{class}} "ash", a list with the following elements(or a  simplified list, if \eqn{onlylogLR=TRUE}, \eqn{minimaloutput=TRUE}   or \eqn{multiseqoutput=TRUE}) \cr
#' \item{fitted.g}{fitted mixture, either a normalmix or unimix}
#' \item{logLR}{log P(D|mle(pi)) - log P(D|null)}
#' \item{loglik}{log P(D|mle(pi))}
#' \item{PosteriorMean}{A vector consisting the posterior mean of beta from the mixture}
#' \item{PosteriorSD}{A vector consisting the corresponding posterior standard deviation}
#' \item{PositiveProb}{A vector of posterior probability that beta is positive}
#' \item{NegativeProb}{A vector of posterior probability that beta is negative}
#' \item{ZeroProb}{A vector of posterior probability that beta is zero}
#' \item{lfsr}{The local false sign rate}
#' \item{lfsra}{The local false sign rate(adjusted)}
#' \item{lfdr}{A vector of estimated local false discovery rate}
#' \item{qvalue}{A vector of q values}
#' \item{fit}{The fitted mixture object by \code{\link{mixEM}} or \code{\link{mixVBEM}} }
#' \item{lambda1}{multiplicative "inflation factor"}
#' \item{lambda2}{additive "inflation factor"}
#' \item{call}{a call in which all of the specified arguments are specified by their full names}
#' \item{data}{a list consisting the input betahat and sebetahat}
#' \item{excludeindex}{the vector of index of observations with 0 standard error; if none, then returns NULL}
#' \item{df}{the specified degrees of freedom for (t) distribution of betahat/sebetahat}
#' \item{model}{either "EE" or "ES", denoting whether exchangeable effects (EE) or exchangeable standardized effects (ES) has been used}
#'
#' @seealso \code{\link{ashci}} for computation of credible intervals after getting the ash object return by \code{ash()}
#' @seealso \code{\link{ashm}} for Multi-model Adaptive Shrinkage function
#'
#' @export
#' @examples 
#' beta = c(rep(0,100),rnorm(100))
#' sebetahat = abs(rnorm(200,0,1))
#' betahat = rnorm(200,beta,sebetahat)
#' beta.ash = ash(betahat, sebetahat)
#' summary(beta.ash)
#' plot(betahat,beta.ash$PosteriorMean,xlim=c(-4,4),ylim=c(-4,4))
#' 
#' CIMatrix=ashci(beta.ash,level=0.95) 
#' print(CIMatrix)
#'
#' betahat=betahat+5
#' beta.ash = ash(betahat, sebetahat)
#' summary(beta.ash)
#' plot(betahat,beta.ash$PosteriorMean)
#'
#' #Testing the non-zero mode feature
#' betahat=betahat+5
#' betan.ash=ash(betahat, sebetahat,nonzeromode=TRUE)
#' plot(betahat, betan.ash$PosteriorMean)
#Things to do:
# check sampling routine
# check number of iterations
ash = function(betahat,sebetahat,method = c("shrink","fdr"), 
               mixcompdist = c("uniform","halfuniform","normal"),
               lambda1=1,lambda2=0,nullcheck=TRUE,df=NULL,randomstart=FALSE,
               nullweight=10,nonzeromode=FALSE, 
               pointmass = FALSE, 
               onlylogLR = FALSE, 
               prior=c("uniform","nullbiased"), 
               mixsd=NULL, VB=FALSE,gridmult=sqrt(2),
               minimaloutput=FALSE,
               multiseqoutput=FALSE,
               g=NULL,
               K=1,
               maxiter = 5000,
               retol=1e-5,
               trace=FALSE,
               cxx=FALSE,
               model=c("EE","ES")
               ){
  
  #method provides a convenient interface to set a particular combinations of parameters for prior an
  #If method is supplied, use it to set up specific values for these parameters; provide warning if values
  #are also specified by user
  #If method is not supplied use the user-supplied values (or defaults if user does not specify them)
  
  squaremK=K #A longer name within the algorithm for easiler debugging purpose.
  if(!missing(method)){
    method = match.arg(method) 
    if(method=="shrink"){
      if(missing(prior)){
        prior = "uniform"
      } else {
        warning("Specification of prior overrides default for method shrink")
      }
      if(missing(pointmass)){
        pointmass=FALSE
      } else {
        warning("Specification of pointmass overrides default for method shrink")
      }
    }
    
    if(method=="fdr"){
      if(missing(prior)){
        prior = "nullbiased"
      } else {
        warning("Specification of prior overrides default for method fdr")
      }
      if(missing(pointmass)){
        pointmass=TRUE
      } else {
        warning("Specification of pointmass overrides default for method fdr")
      }
    }  
  }
  
  if(gridmult<=1&multiseqoutput!=TRUE)
    stop("gridmult must be > 1")
    
  #Dealing with precise input of betahat, currently we exclude them from the EM algorithm
  betahat.input=betahat
  sebetahat.input=sebetahat
  excludeindex=c(1:length(sebetahat.input))[sebetahat.input==0]
  if(length(excludeindex)==0) excludeindex=NULL
  betahat= betahat.input[sebetahat.input!=0]
  sebetahat= sebetahat.input[sebetahat.input!=0]  

  model = match.arg(model)
  if(model=="ES"){ #for ES model, standardize
    betahat = betahat/sebetahat
    sebetahat.orig = sebetahat #store so that can be reinstated later
    sebetahat=1
  }
  
  mixcompdist = match.arg(mixcompdist)
  if(mixcompdist=="normal" & !is.null(df)){
    stop("Error:Normal mixture for student-t likelihood is not yet implemented")
  }
  # if(mixcompdist=="uniform" & pointmass==TRUE){
  #    stop("point mass not yet implemented for uniform or half-uniform")
  #  }
  #  if(mixcompdist=="halfuniform" & pointmass==TRUE){
  #    stop("point mass not yet implemented for uniform or half-uniform")
  #  }
  if(!is.numeric(prior)){
    prior = match.arg(prior)
  }  
  
  if(length(sebetahat)==1){
    sebetahat = rep(sebetahat,length(betahat))
  }
  if(length(sebetahat) != length(betahat)){
    stop("Error: sebetahat must have length 1, or same length as betahat")
  }
  
  completeobs = (!is.na(betahat) & !is.na(sebetahat))
  n=sum(completeobs)
  
  if(n==0){
    if(onlylogLR){
      return(list(pi=NULL, logLR = 0))
    }
    else{
      stop("Error: all input values are missing")
    }
  }
  
  if(missing(trace)){
    if(n>50000){
      trace=TRUE
    }else {trace=FALSE}
  }
  
  if(!is.null(g)){
    maxiter = 1 # if g is specified, don't iterate the EM
    prior = rep(1,ncomp(g)) #prior is not actually used if g specified, but required to make sure EM doesn't produce warning
    null.comp=1 #null.comp also not used, but required 
  } else {
    if(is.null(mixsd)){
      if(nonzeromode){
        mixsd = autoselect.mixsd(betahat[completeobs]-mean(betahat[completeobs]),sebetahat[completeobs],gridmult)
        if(pointmass){ mixsd = c(0,mixsd) }
        nonzeromode.fit=nonzeromodeEM(betahat[completeobs], sebetahat[completeobs], mixsd=mixsd, mixcompdist=mixcompdist,df=df,retol=retol,maxiter=round(maxiter/10),trace=trace,K=squaremK)
        betahat[completeobs]= betahat[completeobs] - nonzeromode.fit$nonzeromode
      }
      else if(nonzeromode & !is.null(df)){
      # stop("Error: Nonzero mean only implemented for df=NULL")
      }
      mixsd = autoselect.mixsd(betahat[completeobs],sebetahat[completeobs],gridmult)
    }
    if(pointmass){
      mixsd = c(0,mixsd)
    }
    
    
    null.comp = which.min(mixsd) #which component is the "null"
    
    k = length(mixsd)
    if(!is.numeric(prior)){
      if(prior=="nullbiased"){ # set up prior to favour "null"
        prior = rep(1,k)
        prior[null.comp] = nullweight #prior 10-1 in favour of null by default
      }else if(prior=="uniform"){
        prior = rep(1,k)
      }
    }
    
    if(length(prior)!=k | !is.numeric(prior)){
      stop("invalid prior specification")
    }
    
    if(randomstart){
      pi = rgamma(k,1,1)
    } else {
      if(k<n){
        pi=rep(1,k)/n #default initialization strongly favours null; puts weight 1/n on everything except null
        pi[null.comp] = (n-k+1)/n #the motivation is data can quickly drive away from null, but tend to drive only slowly toward null.
      } else {
        pi=rep(1,k)/k
      }
    }    
    pi=normalize(pi)
    if(!is.element(mixcompdist,c("normal","uniform","halfuniform"))) stop("Error: invalid type of mixcompdist")
    if(mixcompdist=="normal") g=normalmix(pi,rep(0,k),mixsd)
    if(mixcompdist=="uniform") g=unimix(pi,-mixsd,mixsd)
    if(mixcompdist=="halfuniform"){
      if(min(mixsd)>0){ #simply reflect the components
        g = unimix(c(pi,pi)/2,c(-mixsd,rep(0,k)),c(rep(0,k),mixsd))
        prior = rep(prior, 2)
        pi = rep(pi, 2)
      } else { #define two sets of components, but don't duplicate null component
        null.comp=which.min(mixsd)
        g = unimix(c(pi,pi[-null.comp])/2,c(-mixsd,rep(0,k-1)),c(rep(0,k),mixsd[-null.comp]))
        prior = c(prior,prior[-null.comp])
        pi = c(pi,pi[-null.comp])
      }
    }
  }
  
  pi.fit=EMest(betahat[completeobs],lambda1*sebetahat[completeobs]+lambda2,g,prior,null.comp=null.comp,nullcheck=nullcheck,VB=VB,maxiter = maxiter, cxx=cxx, df=df,trace=trace,K=squaremK)  
  
  #A stringent criteria based on central limit theorem is set to give the user warning message.
  #if(!nonzeromode){
  #  maxsd=max(mixsd)
  #	maxse=quantile(sebetahat[completeobs],0.999)
  #	thresholdval=qnorm(0.999,mean=0,sd=maxse+maxsd)
  #	currentval=abs(sum(betahat[completeobs])/sqrt(length(betahat[completeobs])))
  #	if(currentval>thresholdval){
  #		print("Caution:It's likely that the input data is not coming from a distribution with zero mean, consider to set nonzeromode=TRUE when applying ash()")
  #	}
  #}
  if(!nonzeromode&length(completeobs)>1){
    zvalue=betahat[completeobs]/sebetahat[completeobs]
    abststat=abs(mean(zvalue)/sd(zvalue))*sqrt(length(zvalue))
    if(!is.na(abststat)){
    if(abststat>3.2905){print("Warning: It's likely that the input data is not coming from a distribution with zero mean, consider to set nonzeromode=TRUE when applying ash()")}}
  }
  
  if (!onlylogLR){
      n=length(betahat)
      if (!multiseqoutput){
          ZeroProb = rep(0,length=n)
          NegativeProb = rep(0,length=n)
      }
      if (!minimaloutput){
          PosteriorMean = rep(0,length=n)
          PosteriorSD = rep(0,length=n)
      }
      
            
      if (!multiseqoutput){
          ZeroProb[completeobs] = colSums(comppostprob(pi.fit$g,betahat[completeobs],sebetahat[completeobs],df)[comp_sd(pi.fit$g)==0,,drop=FALSE])
          NegativeProb[completeobs] = cdf_post(pi.fit$g, 0, betahat[completeobs],sebetahat[completeobs],df) - ZeroProb[completeobs]
      }
      if (!minimaloutput){
          PosteriorMean[completeobs] = postmean(pi.fit$g,betahat[completeobs],sebetahat[completeobs],df)
          PosteriorSD[completeobs] = postsd(pi.fit$g,betahat[completeobs],sebetahat[completeobs],df)
      }
      
      
                                        #FOR MISSING OBSERVATIONS, USE THE PRIOR INSTEAD OF THE POSTERIOR
      if (!multiseqoutput){
          ZeroProb[!completeobs] = sum(mixprop(pi.fit$g)[comp_sd(pi.fit$g)==0])
          NegativeProb[!completeobs] = mixcdf(pi.fit$g,0)
          lfsr = compute_lfsr(NegativeProb,ZeroProb)
      }
      if (!minimaloutput){
          PosteriorMean[!completeobs] = calc_mixmean(pi.fit$g)
          PosteriorSD[!completeobs] = calc_mixsd(pi.fit$g)
      }
      if (!minimaloutput & !multiseqoutput){
          PositiveProb = 1- NegativeProb-ZeroProb
          lfsra = compute_lfsra(PositiveProb,NegativeProb,ZeroProb) 
          lfdr = ZeroProb
          qvalue = qval.from.lfdr(lfdr)
      }
  }
  if (!minimaloutput)
      logLR = tail(pi.fit$loglik,1) - pi.fit$null.loglik
  
  if(nonzeromode){
      #Adding back the nonzero mean
      betahat[completeobs]= betahat[completeobs]+nonzeromode.fit$nonzeromode
      if(mixcompdist=="normal"){
        pi.fit$g$mean = rep(nonzeromode.fit$nonzeromode,length(pi.fit$g$pi))
      }
      else if(mixcompdist=="uniform"|mixcompdist=="halfuniform"){
        pi.fit$g$a = pi.fit$g$a + nonzeromode.fit$nonzeromode
        pi.fit$g$b = pi.fit$g$b + nonzeromode.fit$nonzeromode
      }
      PosteriorMean = PosteriorMean + nonzeromode.fit$nonzeromode      
  }	   
  
  if(model=="ES"){
    betahat=betahat*sebetahat.orig
    sebetahat = sebetahat.orig
    PosteriorMean = PosteriorMean * sebetahat
    PosteriorSD= PosteriorSD * sebetahat
  }
  
  loglik = calc_gloglik(pi.fit$g, betahat, sebetahat,df, model) 
  
  if (onlylogLR)
      return(list(fitted.g=pi.fit$g, logLR = logLR, df=df))
  else if (minimaloutput)
      return(list(fitted.g = pi.fit$g, lfsr = lfsr, fit = pi.fit,df=df))
  else if (multiseqoutput)
      return(list(fitted.g = pi.fit$g, logLR = logLR, PosteriorMean = PosteriorMean, PosteriorSD = PosteriorSD, call= match.call(),df=df))
  else{
      result = list(fitted.g = pi.fit$g, logLR = logLR, loglik=loglik, PosteriorMean = PosteriorMean, PosteriorSD = PosteriorSD, PositiveProb = PositiveProb, NegativeProb = NegativeProb, ZeroProb = ZeroProb, lfsr = lfsr,lfsra = lfsra, lfdr = lfdr, qvalue = qvalue, fit = pi.fit, lambda1 = lambda1, lambda2 = lambda2, call = match.call(), data = list(betahat = betahat, sebetahat=sebetahat),excludeindex= excludeindex,df=df,model=model)
      class(result) = "ash"
      return(result)
  }
}

#' @title Function to compute the local false sign rate
#'
#' @param NegativeProb A vector of posterior probability that beta is negative.
#' @param ZeroProb A vector of posterior probability that beta is zero.
#' @return The local false sign rate.
compute_lfsr = function(NegativeProb,ZeroProb){
  ifelse(NegativeProb> 0.5*(1-ZeroProb),1-NegativeProb,NegativeProb+ZeroProb)
}

compute_lfsra = function(PositiveProb, NegativeProb,ZeroProb){
  ifelse(PositiveProb<NegativeProb,2*PositiveProb+ZeroProb,2*NegativeProb+ZeroProb)  
}  

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
#' @param retol the relative tolerance for estimate of nonzero mean,default is 1e-5.
#' @param maxiter the maximum number of iterations performed, default is 500.
#' @param trace a logical variable denoting whether some of the intermediate results of iterations should be displayed to the user, default is FALSE.
#' @param K An integer denoting the order of the SQUAREM scheme. Default is 1,i.e. first-order schemes, which is adequate for most problems. K=2,3 may provide greater speed in some problems, although they are less reliable than the first-order schemes.

#' 
#' @return A list, including the estimates (\eqn{\mu}) and (\eqn{\pi}), the log likelihood for each iteration (NQ)
#' and a flag to indicate convergence.
#' 
#' @export
#' 
#' 
nonzeromodeEM = function(betahat, sebetahat, mixsd, mixcompdist, df=NULL, pi.init=NULL,retol=1e-5,maxiter=500,trace=trace,K=1){
  if(is.null(pi.init)){
    pi.init = rep(1/length(mixsd),length(mixsd))# Use as starting point for pi
  }
  else{
    pi.init=rgamma(length(mixsd),1,1)
  }
  
  if(trace==TRUE){tic()}
  
  #tol in squarem needs special attention as it compares the difference of estimate fixed point (\mu),
  #thus, the absolute tol is not of our interest.here we set it to be relative tol, and our reference 
  #level would be the sample mean, the naive estimator,we keep significant figures to be 6 for fast 
  #convergence failing to set tol to adapt to data would result in keeping the fixed point iteration 
  #running till maxiter, and NaN would be produced, making the method inapplicable
  tol=max(mean(betahat),1)*retol

  if(length(sebetahat)==1){
    sebetahat = rep(sebetahat,length(betahat))
  }
  if(!is.element(mixcompdist,c("normal","uniform","halfuniform"))) stop("Error: invalid type of mixcompdist occcur in nonzeromodeEM()")
  
  if(mixcompdist=="normal" & is.null(df)){
    g=normalmix(pi.init,rep(0,length(mixsd)),mixsd)
    mupi=c(mean(betahat),pi.init)
    res=squarem(par=mupi,fixptfn=nonzeromodeEMfixpoint,objfn=nonzeromodeEMobj,betahat=betahat,sebetahat=sebetahat,mixsd=mixsd,control=list(maxiter=maxiter,tol=tol,trace=trace,K=K))
  }
  else if(mixcompdist=="normal" & !is.null(df)){
    stop("method comp_postsd of normal mixture not yet written for t likelihood")
  }
  else if(mixcompdist=="uniform"){
    g=unimix(pi.init,-mixsd,mixsd)
    mupi=c(mean(betahat),pi.init)    
    res=squarem(par=mupi,fixptfn=nonzeromodeEMoptimfixpoint,objfn=nonzeromodeEMoptimobj,betahat=betahat,sebetahat=sebetahat,g=g,df=df,control=list(maxiter=maxiter,tol=tol,trace=trace,K=K))
  }
  else if(mixcompdist=="halfuniform"){
    g=unimix(c(pi.init, pi.init)/2,c(-mixsd,rep(0,length(mixsd))),c(rep(0,length(mixsd)),mixsd))
    mupi=c(mean(betahat),pi.init/2,pi.init/2)
    res=squarem(par=mupi,fixptfn=nonzeromodeEMoptimfixpoint,objfn=nonzeromodeEMoptimobj,betahat=betahat,sebetahat=sebetahat,g=g,df=df,control=list(maxiter=maxiter,tol=tol,trace=trace,K=K))
  }
  if(trace==TRUE){toc()}
  return(list(nonzeromode=res$par[1],pi=res$par[-1],NQ=-res$value.objfn,niter = res$iter, converged=res$convergence,post=res$par))
}


nonzeromodeEMoptimfixpoint = function(mupi,betahat,sebetahat,g,df){
  mu=mupi[1]
  pimean=normalize(pmax(0,mupi[-1])) #avoid occasional problems with negative pis due to rounding
  matrix_lik=t(compdens_conv(g,betahat-mu,sebetahat,df))
  m=t(pimean * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum=rowSums(m)
  classprob=m/m.rowsum #an n by k matrix
  pinew=normalize(colSums(classprob))
  munew=optimize(f=nonzeromodeEMoptim,interval=c(min(betahat),max(betahat)), pinew=pinew,betahat=betahat,sebetahat=sebetahat,g=g,df=df)$minimum
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
#' @param matrix_lik a n by k matrix with (j,k)th element equal to \eqn{f_k(x_j)}.
#' @param prior a k vector of the parameters of the Dirichlet prior on \eqn{\pi}. Recommended to be rep(1,k)
#' @param pi.init the initial value of the posterior parameters. If not specified defaults to the prior parameters.
#' @param tol the tolerance for convergence of log-likelihood bound.
#' @param maxiter the maximum number of iterations performed
#' @param trace a logical variable denoting whether some of the intermediate results of iterations should be displayed to the user. Default is FALSE.
#' @param K An integer denoting the order of the SQUAREM scheme. Default is 1,i.e. first-order schemes, which is adequate for most problems. K=2,3 may provide greater speed in some problems, although they are less reliable than the first-order schemes.
#' 
#' @return A list, whose components include point estimates (pihat), 
#' the parameters of the fitted posterior on \eqn{\pi} (pipost),
#' the bound on the log likelihood for each iteration (B)
#' and a flag to indicate convergence (converged).
#'  
#' @export
#' 
mixVBEM = function(matrix_lik, prior, pi.init = NULL,tol=1e-7, maxiter=5000,trace=FALSE,K=1){
  k=ncol(matrix_lik)
  if(is.null(pi.init)){
    pi.init = rep(1,k)# Use as starting point for pi
  } 
  res = squarem(par=pi.init,fixptfn=VBfixpoint, objfn=VBnegpenloglik,matrix_lik=matrix_lik, prior=prior, control=list(maxiter=maxiter,tol=tol,trace=trace,K=K))
  return(list(pihat = res$par/sum(res$par), B=res$value.objfn, niter = res$iter, converged=res$convergence,post=res$par))
}


VBfixpoint = function(pipost, matrix_lik, prior){  
  n=nrow(matrix_lik)
  k=ncol(matrix_lik)
  avgpipost = matrix(exp(rep(digamma(pipost),n)-rep(digamma(sum(pipost)),k*n)),ncol=k,byrow=TRUE)
  classprob = avgpipost*matrix_lik
  classprob = classprob/rowSums(classprob) # n by k matrix
  pipostnew = colSums(classprob) + prior
  return(pipostnew)
}

VBnegpenloglik=function(pipost,matrix_lik,prior){
  return(-VBpenloglik(pipost,matrix_lik,prior))
}

VBpenloglik = function(pipost, matrix_lik, prior){
  n=nrow(matrix_lik)
  k=ncol(matrix_lik)
  avgpipost = matrix(exp(rep(digamma(pipost),n)-rep(digamma(sum(pipost)),k*n)),ncol=k,byrow=TRUE)
  classprob = avgpipost*matrix_lik
  classprob = classprob/rowSums(classprob) # n by k matrix
  
  B= sum(classprob*log(avgpipost*matrix_lik),na.rm=TRUE) - diriKL(prior,pipost) - sum(classprob*log(classprob)) 
  return(B)
}



#' @title Estimate mixture proportions of a mixture model by EM algorithm
#'
#' @description Given the individual component likelihoods for a mixture model, estimates the mixture proportions by an EM algorithm.
#'
#' @details Fits a k component mixture model \deqn{f(x|\pi)= \sum_k \pi_k f_k(x)} to independent
#' and identically distributed data \eqn{x_1,\dots,x_n}. 
#' Estimates mixture proportions \eqn{\pi} by maximum likelihood, or by maximum a posteriori (MAP) estimation for a Dirichlet prior on \eqn{\pi} 
#' (if a prior is specified).  Uses the SQUAREM package to accelerate convergence of EM. Used by the ash main function; there is no need for a user to call this 
#' function separately, but it is exported for convenience.
#'
#' 
#' @param matrix_lik, a n by k matrix with (j,k)th element equal to \eqn{f_k(x_j)}.
#' @param prior, a k vector of the parameters of the Dirichlet prior on \eqn{\pi}. Recommended to be rep(1,k)
#' @param pi.init, the initial value of \eqn{\pi} to use. If not specified defaults to (1/k,...,1/k).
#' @param tol, the tolerance for convergence of log-likelihood.
#' @param maxiter the maximum number of iterations performed
#' @param trace a logical variable denoting whether some of the intermediate results of iterations should be displayed to the user. Default is FALSE.
#' @param K An integer denoting the order of the SQUAREM scheme. Default is 1,i.e. first-order schemes, which is adequate for most problems. K=2,3 may provide greater speed in some problems, although they are less reliable than the first-order schemes.
#' 
#' @return A list, including the estimates (pihat), the log likelihood for each interation (B)
#' and a flag to indicate convergence
#'  
#' @export
#' 
#' 
mixEM = function(matrix_lik, prior, pi.init = NULL,tol=1e-7, maxiter=5000,trace=FALSE,K=1){
  k=dim(matrix_lik)[2]
  if(is.null(pi.init)){
    pi.init = rep(1/k,k)# Use as starting point for pi
  } 
  res = squarem(par=pi.init,fixptfn=fixpoint, objfn=negpenloglik,matrix_lik=matrix_lik, prior=prior, control=list(maxiter=maxiter,tol=tol,trace=trace,K=K))
  return(list(pihat = normalize(pmax(0,res$par)), B=res$value.objfn, 
              niter = res$iter, converged=res$convergence))
}

# helper functions used by mixEM
normalize = function(x){return(x/sum(x))}

fixpoint = function(pi, matrix_lik, prior){  
  pi = normalize(pmax(0,pi)) #avoid occasional problems with negative pis due to rounding
  m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  classprob = m/m.rowsum #an n by k matrix
  pinew = normalize(colSums(classprob) + prior - 1)
  return(pinew)
}

negpenloglik = function(pi,matrix_lik,prior){return(-penloglik(pi,matrix_lik,prior))}

penloglik = function(pi, matrix_lik, prior){
  pi = normalize(pmax(0,pi))
  m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  loglik = sum(log(m.rowsum))
  subset = (prior != 1.0)
  priordens = sum((prior-1)[subset]*log(pi[subset]))
  return(loglik+priordens)
}

#The kth element of this vector is the derivative 
#of the loglik for $\pi=(\pi_0,...,1-\pi_0,...)$ with respect to $\pi_0$ at $\pi_0=1$.
gradient = function(matrix_lik){
  n = nrow(matrix_lik)
  grad = n - colSums(matrix_lik/matrix_lik[,1]) 
  return(grad)
}

# mixEM = function(matrix_lik, prior, pi.init = NULL,tol=0.0001, maxiter=5000){
#   n=nrow(matrix_lik)
#   k=ncol(matrix_lik)
#   B = rep(0,maxiter)
#   pi = pi.init
#   if(is.null(pi.init)){
#     pi = rep(1/k,k)# Use as starting point for pi
#   } 
#   pi = ifelse(pi<1e-5,1e-5,pi) #set any estimates that are too small to be just very small
#   pi = normalize(pi)
#   
#   loglik = rep(0,maxiter)
#   priordens= rep(0,maxiter)
#   m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
#   m.rowsum = rowSums(m)
#   loglik[1] = sum(log(m.rowsum))
#   priordens[1] = sum((prior-1)*log(pi)) 
#   classprob = m/m.rowsum #an n by k matrix
#   i=1
#   if(maxiter >= 2){
#     for(i in 2:maxiter){  
#       pi = colSums(classprob) + prior-1
#       pi = ifelse(pi<1e-5,1e-5,pi) #set any estimates that are less than zero, which can happen with prior<1, to 0
#       pi = normalize(pi)
#         
#       #Now re-estimate pi
#       m  = t(pi * t(matrix_lik)) 
#       m.rowsum = rowSums(m)
#       loglik[i] = sum(log(m.rowsum))
#       priordens[i] = sum((prior-1)*log(pi)) 
#       classprob = m/m.rowsum
#     
#     
#       if(abs(loglik[i]+priordens[i]-loglik[i-1]-priordens[i-1])<tol) break;
#     }
#   }
#   converged=(abs(loglik[i]+priordens[i]-loglik[i-1]-priordens[i-1])<tol)
#   if(!converged){
#       warning("EM algorithm in function mixEM failed to converge. Results may be unreliable. Try increasing maxiter and rerunning.")
#   }
#   return(list(pihat = pi, B=loglik[1:i], 
#               niter = i, converged=converged))
# }


#' @title estimate mixture proportions of sigmaa by EM algorithm
#'
#' @param betahat (n vector of observations) 
#' @param sebetahat (n vector of standard errors/deviations of observations)
#' @param g the prior distribution for beta (usually estimated from the data; this is used primarily in simulated data to do computations with the "true" g)
#' @param prior string, or numeric vector indicating Dirichlet prior on mixture proportions (defaults to "uniform", or (1,1...,1); also can be "nullbiased" (nullweight,1,...,1) to put more weight on first component)
#' @param null.comp the position of the null component
#' @param nullcheck  whether to check that any fitted model exceeds the "null" likelihood
#' in which all weight is on the first component
#' @param VB whether to use Variational Bayes to estimate mixture proportions (instead of EM to find MAP estimate), see \code{\link{mixVBEM}} and \code{\link{mixEM}}
#' @param maxiter maximum number of iterations of the EM algorithm.
#' @param cxx flag to indicate whether to use the c++ (Rcpp) version. After application of Squared extrapolation methods for accelerating fixed-point iterations (R Package "SQUAREM"), the c++ version is no longer faster than non-c++ version, thus we do not recommend using this one, and might be removed at any point. 
#' @param df appropriate degrees of freedom for (t) distribution of betahat/sebetahat, default is NULL(Gaussian)
#' @param trace a logical variable denoting whether some of the intermediate results of iterations should be displayed to the user. Default is FALSE.
#' @param K An integer denoting the order of the SQUAREM scheme. Default is 1,i.e. first-order schemes, which is adequate for most problems. K=2,3 may provide greater speed in some problems, although they are less reliable than the first-order schemes.
#' @return A list, including the final loglikelihood, the null loglikelihood, a n by k likelihoodmatrix with (j,k)th element equal to \eqn{f_k(x_j)},and a flag to indicate convergence.
#
#prior gives the parameter of a Dirichlet prior on pi
#(prior is used to encourage results towards smallest value of sigma when
#likelihood is flat)
#nullcheck indicates whether to check whether the loglike exceeds the null
#(may not want to use if prior is used)
#VB provides an approach to estimate the approximate posterior distribution
#of mixture proportions of sigmaa by variational Bayes method
#(use Dirichlet prior and approximate Dirichlet posterior)
#if cxx TRUE use cpp version of R function mixEM

EMest = function(betahat,sebetahat,g,prior,null.comp=1,nullcheck=TRUE,VB=FALSE, maxiter=5000, cxx=FALSE, df=NULL,trace=FALSE,K=1){ 
  
  
  pi.init = g$pi
  k=ncomp(g)
  n = length(betahat)
  tol = min(0.1/n,1e-5) # set convergence criteria to be more stringent for larger samples
  
  if(trace==TRUE){tic()}

  matrix_lik = t(compdens_conv(g,betahat,sebetahat,df))

  
  #checks whether the gradient at pi0=1 is positive (suggesting that this is a fixed point)
  #if(nullcheck){
  #  if(all(gradient(matrix_lik)>=0)){
  #    pi.init=rep(0,k)
  #    pi.init[null.comp]=1 #this will make pi.init=(1,0,0...,0) which is a fixed point of the EM
  #  }
  #}
  
  if(VB==TRUE){
    EMfit=mixVBEM(matrix_lik,prior,maxiter=maxiter,trace,K=K)}
  else{
    if (cxx==TRUE){
      EMfit = cxxMixEM(matrix_lik,prior,pi.init,1e-5, maxiter) #currently use different convergence criteria for cxx version 
      if(!EMfit$converged){
        warning("EM algorithm in function cxxMixEM failed to converge. Results may be unreliable. Try increasing maxiter and rerunning.")
      }
    }
    else{
      EMfit = mixEM(matrix_lik,prior,pi.init,tol, maxiter,trace,K=K)
      if(!EMfit$converged & !(maxiter==1)){
        warning("EM algorithm in function mixEM failed to converge. Results may be unreliable. Try increasing maxiter and rerunning.")
      }
    }
  }
  
  pi = EMfit$pihat     
  penloglik = EMfit$B 
  converged = EMfit$converged
  niter = EMfit$niter
  
  loglik.final =  penloglik(pi,matrix_lik,1) #compute penloglik without penalty
  null.loglik = sum(log(matrix_lik[,null.comp]))  
  
  if(nullcheck==TRUE & VB==FALSE){ #null check doesn't work with VB yet
    pinull = rep(0,k)
    pinull[null.comp]=1
    null.penloglik = penloglik(pinull,matrix_lik,prior)
    final.penloglik = penloglik(pi,matrix_lik,prior)
    
    if(null.penloglik > final.penloglik){ #check whether exceeded "null" likelihood where everything is null
      pi=pinull
      loglik.final=penloglik(pi,matrix_lik,1)
    }
  }
  
  g$pi=pi
  if(trace==TRUE){toc()}
  
  return(list(loglik=loglik.final,null.loglik=null.loglik,
              matrix_lik=matrix_lik,converged=converged,g=g))
}




#' @title Compute Posterior
#'
#' @description Return the posterior on beta given a prior (g) that is a mixture of normals (class normalmix) 
#' and observation \eqn{betahat ~ N(beta,sebetahat)}
#'
#' @details This can be used to obt
#'
#' @param g a normalmix with components indicating the prior; works only if g has means 0
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
  sebetahat=sebetahat[sebetahat!=0] #To avoid exact measure causing (usually by mistake)
  sigmaamin = min(sebetahat)/10 #so that the minimum is small compared with measurement precision
  if(all(betahat^2<=sebetahat^2)){
    sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
  }else{
    sigmaamax = 2*sqrt(max(betahat^2-sebetahat^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2   
  }
  if(mult==0){
    return(c(0,sigmaamax/2))
  }else{
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
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
VB.update = function(matrix_lik,pipost,prior){
  n=dim(matrix_lik)[1]
  k=dim(matrix_lik)[2]
  avgpipost = matrix(exp(rep(digamma(pipost),n)-rep(digamma(sum(pipost)),k*n)),ncol=k,byrow=TRUE)
  classprob = avgpipost * matrix_lik
  classprob = classprob/rowSums(classprob) # n by k matrix
  B = sum(classprob*log(avgpipost*matrix_lik),na.rm=TRUE) - diriKL(prior,pipost) #negative free energy
  return(list(classprob=classprob,B=B))
}

#Helper Function for nonzeromodeEM, from MATLAB Package
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
   type <- match.arg(type)
   assign(".type", type, envir=baseenv())
   if(gcFirst) gc(FALSE)
   tic <- proc.time()[type]         
   assign(".tic", tic, envir=baseenv())
   invisible(tic)
}

toc <- function()
{
   type <- get(".type", envir=baseenv())
   toc <- proc.time()[type]
   tic <- get(".tic", envir=baseenv())
   print(toc - tic)
   invisible(toc)
}




#Not Used #' @title Faster version of function ash
# #'
# #' @description This function has similar functionality as ash, but only returns some of the outputs.
# #'
# #' @param betahat, a p vector of estimates
# #' @param sebetahat, a p vector of corresponding standard errors
# #' @param nullcheck: whether to check that any fitted model exceeds the "null" likelihood in which all weight is on the first component
# #' @param randomstart: logical, indicating whether to initialize EM randomly. If FALSE, then initializes to prior mean (for EM algorithm) or prior (for VBEM)
# #' @param pointmass: logical, indicating whether to use a point mass at zero as one of components for a mixture distribution
# #' @param onlylogLR: logical, indicating whether to use this function to get logLR. Skip posterior prob, posterior mean, lfdr...
# #' @param prior: string, or numeric vector indicating Dirichlet prior on mixture proportions (defaults to "uniform", or 1,1...,1; also can be "nullbiased" 1,1/k-1,...,1/k-1 to put more weight on first component)
# #' @param mixsd: vector of sds for underlying mixture components
# #' @param VB: whether to use Variational Bayes to estimate mixture proportions (instead of EM to find MAP estimate)
# #' @param gridmult: the multiplier by which the default grid values for mixsd differ by one another. (Smaller values produce finer grids)
# #' @param g: the prior distribution for beta (usually estimated from the data; this is used primarily in simulated data to do computations with the "true" g)
# #' @param cxx: flag to indicate whether to use the c++ (Rcpp) version
# #'
# #' @return a list with elements fitted.g is fitted mixture
# #' logLR : logP(D|mle(pi)) - logP(D|null)
# #'
# #' @export
# fast.ash = function(betahat,sebetahat, 
#                     nullcheck=TRUE,randomstart=FALSE, 
#                     pointmass = TRUE,    
#                     prior=c("nullbiased","uniform"), 
#                     mixsd=NULL, VB=FALSE,gridmult=4,
#                     g=NULL, cxx=TRUE,
#                     onlylogLR = FALSE,df=NULL){
#   
#   if(onlylogLR){
#     pointmass <- TRUE  
#   }
#   
#   #If method is supplied, use it to set up defaults; provide warning if these default values
#   #are also specified by user
#   if(!is.numeric(prior)){
#     prior = match.arg(prior)
#   }
#   
#   if(length(sebetahat)==1){
#     sebetahat = rep(sebetahat,length(betahat))
#   }
#   if(length(sebetahat) != length(betahat)){
#     stop("Error: sebetahat must have length 1, or same length as betahat")
#   }
#   
#   completeobs = (!is.na(betahat) & !is.na(sebetahat))
#   if(sum(completeobs)==0){
#     if(onlylogLR){
#       return(list(pi=NULL, logLR = 0))
#     }else{
#       stop("Error: all input values are missing")
#     }
#   }  
#   
#   if(is.null(mixsd)){
#     mixsd= autoselect.mixsd(betahat[completeobs],sebetahat[completeobs],gridmult)
#   }
#   if(pointmass){
#     mixsd = c(0,mixsd)
#   }
#   
#   k=length(mixsd)  
#   null.comp = which.min(mixsd) #which component is the "null"
#   
#   if(!is.numeric(prior)){
#     if(prior=="nullbiased"){ # set up prior to favour "null"
#       prior = rep(1,k)
#       prior[null.comp] = 10 #prior 10-1 in favour of null
#     }else if(prior=="uniform"){
#       prior = rep(1,k)
#     }
#   }
#   
#   if(length(prior)!=k | !is.numeric(prior)){
#     stop("invalid prior specification")
#   }
#   
#   if(missing(g)){
#     pi = prior^2 #default is to initialize pi at prior (mean)
#     if(randomstart){pi=rgamma(k,1,1)}
#     pi=normalize(pi)
#     g=normalmix(pi,rep(0,k),mixsd)
#     maxiter = 5000
#   } else {
#     maxiter = 1; # if g is specified, don't iterate the EM 
#   }
#   
#   pi.fit=EMest(betahat[completeobs],sebetahat[completeobs],g,prior,null.comp=null.comp,nullcheck=nullcheck,VB=VB,maxiter = maxiter, cxx=cxx, df=df)  
#   
#   if(onlylogLR){
#     logLR = tail(pi.fit$loglik,1) - pi.fit$null.loglik
#     return(list(pi=pi.fit$pi, logLR = logLR))
#   }else{
#     
#     n=length(betahat)
#     PosteriorMean = rep(0,length=n)
#     PosteriorSD=rep(0,length=n)
#     
#     if(is.null(df)){
#       PosteriorMean[completeobs] = postmean(pi.fit$g,betahat[completeobs],sebetahat[completeobs])
#       PosteriorSD[completeobs] =postsd(pi.fit$g,betahat[completeobs],sebetahat[completeobs]) 
#     }
#     else{
#       PosteriorMean[completeobs] = postmean_t(pi.fit$g,betahat[completeobs],sebetahat[completeobs],df)
#       PosteriorSD[completeobs] =postsd_t(pi.fit$g,betahat[completeobs],sebetahat[completeobs],df)
#     }
#     #FOR MISSING OBSERVATIONS, USE THE PRIOR INSTEAD OF THE POSTERIOR
#     PosteriorMean[!completeobs] = mixmean(pi.fit$g)
#     PosteriorSD[!completeobs] =mixsd(pi.fit$g)  
#     
#     result = list(fitted.g=pi.fit$g,PosteriorMean = PosteriorMean,PosteriorSD=PosteriorSD,call=match.call(),data=list(betahat = betahat, sebetahat=sebetahat))
#     return(result)
#   }
#   #if(nsamp>0){
#   #  sample = posterior_sample(post,nsamp)
#   #}
# }
