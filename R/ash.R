#' @useDynLib ashr
#' @import truncnorm SQUAREM doParallel pscl Rcpp
#
#


#' @title Main Adaptive Shrinkage function
#'
#' @description Takes vectors of estimates (betahat) and their standard errors (sebetahat), and applies
#' shrinkage to them, using Empirical Bayes methods, to compute shrunk estimates for beta.
#'
#' @details This function is actually just a simple wrapper that passes its parameters to ash.workhorse which provides more documented options for advanced use. See readme for more details. 
#' 
#' @param betahat  a p vector of estimates 
#' @param sebetahat a p vector of corresponding standard errors
#' @param mixcompdist distribution of components in mixture ( "uniform","halfuniform" or "normal"), the default is "uniform". If you believe your effects may be asymmetric, use "halfuniform". The use of "normal" is permitted only if df=NULL.
#' @param df appropriate degrees of freedom for (t) distribution of betahat/sebetahat, default is NULL which is actually treated as infinity (Gaussian)
#' 
#' @return ash returns an object of \code{\link[base]{class}} "ash", a list with some or all of the following elements (determined by outputlevel) \cr
#' \item{fitted.g}{fitted mixture, either a normalmix or unimix}
#' \item{loglik}{log P(D|mle(pi))}
#' \item{PosteriorMean}{A vector consisting the posterior mean of beta from the mixture}
#' \item{PosteriorSD}{A vector consisting the corresponding posterior standard deviation}
#' \item{PositiveProb}{A vector of posterior probability that beta is positive}
#' \item{NegativeProb}{A vector of posterior probability that beta is negative}
#' \item{ZeroProb}{A vector of posterior probability that beta is zero}
#' \item{lfsr}{The local false sign rate}
#' \item{lfdr}{A vector of estimated local false discovery rate}
#' \item{qvalue}{A vector of q values}
#' \item{fit}{The fitted mixture object}
#' \item{call}{a call in which all of the specified arguments are specified by their full names}
#' \item{excludeindex}{the vector of index of observations with 0 standard error; if none, then returns NULL}
#' \item{model}{either "EE" or "ET", denoting whether exchangeable effects (EE) or exchangeable T stats (ET) has been used}
#' \item{optmethod}{the optimization method used}
#' \item{data}{a list consisting the input betahat and sebetahat (only included if outputlevel>2)}
#' 
#' @seealso \code{\link{ash.workhorse}} for complete specification of ash function
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
#' #Illustrating the non-zero mode feature
#' betahat=betahat+5
#' beta.ash = ash(betahat, sebetahat)
#' plot(betahat,beta.ash$PosteriorMean)
#' summary(beta.ash)
#' betan.ash=ash(betahat, sebetahat,nonzeromode=TRUE)
#' plot(betahat, betan.ash$PosteriorMean)
#' summary(betan.ash)
ash = function(betahat,sebetahat,mixcompdist = c("uniform","halfuniform","normal"),df=NULL,...){
  return(modifyList(ash.workhorse(betahat,sebetahat,mixcompdist=mixcompdist,df=df,...),list(call=match.call())))
}


#' @title Detailed Adaptive Shrinkage function
#'
#' @description Takes vectors of estimates (betahat) and their standard errors (sebetahat), and applies
#' shrinkage to them, using Empirical Bayes methods, to compute shrunk estimates for beta. This is the more detailed version of ash for "research" use. 
#' Most users will be happy with the ash function, which provides the same usage, but documents only the main options for simplicity. 
#'
#' @details See readme for more details.
#' 
#' @param betahat  a p vector of estimates 
#' @param sebetahat a p vector of corresponding standard errors
#' @param method specifies how ash is to be run. Can be "shrinkage" (if main aim is shrinkage) or "fdr" (if main aim is to assess fdr or fsr)
#' This is simply a convenient way to specify certain combinations of parameters: "shrinkage" sets pointmass=FALSE and prior="uniform";
#' "fdr" sets pointmass=TRUE and prior="nullbiased".
#' @param mixcompdist distribution of components in mixture ( "uniform","halfuniform" or "normal"), the default value is "uniform"
#'
#' @param optmethod specifies optimization method used. Default is "mixIP", an interior point method, if REBayes is installed; otherwise a slower EM algorithm is used.
#' @param df appropriate degrees of freedom for (t) distribution of betahat/sebetahat, default is NULL(Gaussian)
#' @param nullweight scalar, the weight put on the prior under "nullbiased" specification, see \code{prior}
#' @param randomstart logical, indicating whether to initialize EM randomly. If FALSE, then initializes to prior mean (for EM algorithm) or prior (for VBEM)
#' @param nonzeromode logical, indicating whether to use a non-zero unimodal mixture(default is "FALSE")
#' @param pointmass logical, indicating whether to use a point mass at zero as one of components for a mixture distribution
#' @param prior string, or numeric vector indicating Dirichlet prior on mixture proportions (defaults to "uniform", or (1,1...,1); also can be "nullbiased" (nullweight,1,...,1) to put more weight on first component)
#' @param mixsd vector of sds for underlying mixture components 
#' @param gridmult the multiplier by which the default grid values for mixsd differ by one another. (Smaller values produce finer grids)
#' @param outputlevel determines amount of output [0=just fitted g; 1=also PosteriorMean and PosteriorSD; 2= all but data; 3=everything, including copy of input data]
#' @param g the prior distribution for beta (usually estimated from the data; this is used primarily in simulated data to do computations with the "true" g)
#' @param fixg if TRUE, don't estimate g but use the specified g - useful for computations under the "true" g in simulations
#' @param VB (deprecated, use optmethod) whether to use Variational Bayes to estimate mixture proportions (instead of EM to find MAP estimate), see \code{\link{mixVBEM}} and \code{\link{mixEM}}
#' @param cxx flag (deprecated, use optmethod) to indicate whether to use the c++ (Rcpp) version. After application of Squared extrapolation methods for accelerating fixed-point iterations (R Package "SQUAREM"), the c++ version is no longer faster than non-c++ version, thus we do not recommend using this one, and might be removed at any point. 
#' @param model c("EE","ET") specifies whether to assume exchangeable effects (EE) or exchangeable T stats (ET).
#' @param control A list of control parameters for the optmization algorithm. Default value is set to be   control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE). User may supply changes to this list of parameter, say, control=list(maxiter=10000,trace=TRUE)

#' 
#' @return ash returns an object of \code{\link[base]{class}} "ash", a list with some or all of the following elements (determined by outputlevel) \cr
#' \item{fitted.g}{fitted mixture, either a normalmix or unimix}
#' \item{loglik}{log P(D|mle(pi))}
#' \item{PosteriorMean}{A vector consisting the posterior mean of beta from the mixture}
#' \item{PosteriorSD}{A vector consisting the corresponding posterior standard deviation}
#' \item{PositiveProb}{A vector of posterior probability that beta is positive}
#' \item{NegativeProb}{A vector of posterior probability that beta is negative}
#' \item{ZeroProb}{A vector of posterior probability that beta is zero}
#' \item{lfsr}{The local false sign rate}
#' \item{lfdr}{A vector of estimated local false discovery rate}
#' \item{qvalue}{A vector of q values}
#' \item{fit}{The fitted mixture object}
#' \item{call}{a call in which all of the specified arguments are specified by their full names}
#' \item{excludeindex}{the vector of index of observations with 0 standard error; if none, then returns NULL}
#' \item{model}{either "EE" or "ET", denoting whether exchangeable effects (EE) or exchangeable T stats (ET) has been used}
#' \item{optmethod}{the optimization method used}
#' \item{data}{a list consisting the input betahat and sebetahat (only included if outputlevel>2)}
#'
#' @seealso \code{\link{ash}} for simplified specification of ash function
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
#' #Testing the non-zero mode feature
#' betahat=betahat+5
#' beta.ash = ash(betahat, sebetahat)
#' plot(betahat,beta.ash$PosteriorMean)
#' summary(beta.ash)
#' betan.ash=ash(betahat, sebetahat,nonzeromode=TRUE)
#' plot(betahat, betan.ash$PosteriorMean)
#' summary(betan.ash)
#' 
#' #Running ash with a pre-specified g, rather than estimating it
#' beta = c(rep(0,100),rnorm(100))
#' sebetahat = abs(rnorm(200,0,1))
#' betahat = rnorm(200,beta,sebetahat)
#' true_g = normalmix(c(0.5,0.5),c(0,0),c(0,1)) # define true g 
#' #Passing this g into ash causes it to i) take the sd and the means for each component from this g, and ii) initialize pi to the value from this g.
#' beta.ash = ash(betahat, sebetahat,g=true_g,fixg=TRUE)
ash.workhorse = function(betahat,sebetahat,
               method = c("fdr","shrink"),
               mixcompdist = c("uniform","halfuniform","normal"),
               optmethod = c("mixIP","cxxMixSquarem","mixEM","mixVBEM"),
               df=NULL,randomstart=FALSE,
               nullweight=10,nonzeromode=FALSE, 
               pointmass = TRUE, 
               prior=c("nullbiased","uniform"), 
               mixsd=NULL, gridmult=sqrt(2),
               outputlevel=2,
               g=NULL,
               fixg=FALSE,
               cxx=FALSE,
               VB=FALSE,
               model=c("EE","ET"),
               control=list()
){
  
  ##1.Handling Input Parameters
  

  # Set optimization method (optmethod)
  
  # if user tries to set both optmethod and VB/cxx that's an error
  if(!missing(optmethod) && (!missing(VB) || !missing(cxx)) ){
    stop("VB and cxx options are deprecated and incompatible with optmethod; use optmethod instead")
  }

  # if no optmethod specified, select a default
  if(missing(optmethod)){
    if(require(REBayes,quietly=TRUE)){ #check whether REBayes package is present
      optmethod = "mixIP"
    } else{  #If REBayes package missing
      message("Due to absence of package REBayes, switching to EM algorithm")
      if(require(Rcpp)){
        optmethod = "cxxMixSquarem"}
      else {
        optmethod = "mixEM" #fallback if neither Rcpp or REBayes are installed
        message("Using vanilla EM; for faster performance install REBayes (preferred) or Rcpp")
      }
    }
  } else { #if optmethod specified
    optmethod = match.arg(optmethod)
  }
  
  if(!missing(VB)){
    warning("VB option is deprecated, use optmethod instead")
    if(VB==TRUE){optmethod = "mixVBEM"}
  }
    
  if(!missing(cxx)){
    warning("cxx option is deprecated, use optmethod instead")
    if (cxx==TRUE){optmethod = "cxxMixSquarem"}
    if (cxx==FALSE){optmethod = "mixEM"}
  }
    
  if(optmethod == "mixIP"){assertthat::assert_that(require(REBayes,quietly=TRUE))}
  if(optmethod == "cxxMixSquarem"){assertthat::assert_that(require(Rcpp))}
  
  #method provides a convenient interface to set a particular combinations of parameters for prior an
  #If method is supplied, use it to set up specific values for these parameters; provide warning if values
  #are also specified by user
  #If method is not supplied use the user-supplied values (or defaults if user does not specify them)  
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
  
  #Dealing with precise input of betahat, currently we exclude them from the EM algorithm
  betahat.input=betahat
  sebetahat.input=sebetahat
  excludeindex=c(1:length(sebetahat.input))[sebetahat.input==0]
  if(length(excludeindex)==0) excludeindex=NULL
  betahat= betahat.input[sebetahat.input!=0]
  sebetahat= sebetahat.input[sebetahat.input!=0]
  
  #Set observations with infinite standard errors to missing
  #later these missing observations will be ignored in EM, and posterior will be same as prior.
  sebetahat[sebetahat==Inf]=NA 
  betahat[sebetahat==Inf]=NA
  
  model = match.arg(model)
  if(model=="ET"){ #for ET model, standardize
    betahat = betahat/sebetahat
    sebetahat.orig = sebetahat #store so that can be reinstated later
    sebetahat=1
  }
  
  mixcompdist = match.arg(mixcompdist)
  if(!is.numeric(prior)){  prior = match.arg(prior)  } 
  
  if(mixcompdist=="normal" & !is.null(df)){
    stop("Error:Normal mixture for student-t likelihood is not yet implemented")
  }
  if(mixcompdist=="halfuniform" & prior!="nullbiased"){
    warning("Use of halfuniform without nullbiased prior can lead to misleading local false sign rates, and so is not recommended")
  }
  
  if(gridmult<=1){stop("gridmult must be > 1")}  
  if(length(sebetahat)==1){  sebetahat = rep(sebetahat,length(betahat))  }
  if(length(sebetahat) != length(betahat)){
    stop("Error: sebetahat must have length 1, or same length as betahat")
  }
  
  completeobs = (!is.na(betahat) & !is.na(sebetahat))
  n=sum(completeobs)
  
  #Handling control variables
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  if(n>50000){control.default$trace=TRUE}
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)
  if(controlinput$maxiter==0){
    stop("option control$maxiter=0 deprecated; used fixg=TRUE instead")
  }
  
  if(n==0){
    stop("Error: all input values are missing") 
  }
  
  
  ##2. Generating mixture distribution
  
  if(fixg & missing(g)){stop("if fixg=TRUE then you must specify g!")}
  
  if(!is.null(g)){
    k=ncomp(g)
    null.comp=1 #null.comp not actually used unless randomstart true 
    prior = setprior(prior,k,nullweight,null.comp)
    if(randomstart){pi = initpi(k,n,null.comp,randomstart)
                    g$pi=pi} #if g specified, only initialize pi if randomstart is TRUE 
  } else {
    if(is.null(mixsd)){
      if(nonzeromode){
        mixsd = autoselect.mixsd(betahat[completeobs]-mean(betahat[completeobs]),sebetahat[completeobs],gridmult)
        if(pointmass){ mixsd = c(0,mixsd) }
        nonzeromode.fit=nonzeromodeEM(betahat[completeobs], sebetahat[completeobs], mixsd=mixsd, mixcompdist=mixcompdist,df=df,control= controlinput)
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
    prior = setprior(prior,k,nullweight,null.comp)
    pi = initpi(k,n,null.comp,randomstart)
 
    
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
  
  
  ##3. Fitting the mixture
  if(!fixg){
    pi.fit=estimate_mixprop(betahat[completeobs],sebetahat[completeobs],g,prior,null.comp=null.comp,
               optmethod=optmethod,df=df,control=controlinput)  
  } else {
    pi.fit = list(g=g)
  }
  
  ##4. Computing the posterior
  
  n = length(betahat)
  
  if (outputlevel > 0) {
    PosteriorMean = rep(0,length = n)
    PosteriorSD = rep(0,length = n)
    PosteriorMean[completeobs] = postmean(pi.fit$g,betahat[completeobs],sebetahat[completeobs],df)
    PosteriorSD[completeobs] = postsd(pi.fit$g,betahat[completeobs],sebetahat[completeobs],df)
    #FOR MISSING OBSERVATIONS, USE THE PRIOR INSTEAD OF THE POSTERIOR
    PosteriorMean[!completeobs] = calc_mixmean(pi.fit$g)
    PosteriorSD[!completeobs] = calc_mixsd(pi.fit$g)
  }
  if (outputlevel > 1) {
    ZeroProb = rep(0,length = n)
    NegativeProb = rep(0,length = n)
    ZeroProb[completeobs] = colSums(comppostprob(pi.fit$g,betahat[completeobs],sebetahat[completeobs],df)[comp_sd(pi.fit$g) ==
                                                                                                            0,,drop = FALSE])
    NegativeProb[completeobs] = cdf_post(pi.fit$g, 0, betahat[completeobs],sebetahat[completeobs],df) - ZeroProb[completeobs]
    #FOR MISSING OBSERVATIONS, USE THE PRIOR INSTEAD OF THE POSTERIOR
    ZeroProb[!completeobs] = sum(mixprop(pi.fit$g)[comp_sd(pi.fit$g) == 0])
    NegativeProb[!completeobs] = mixcdf(pi.fit$g,0)
    lfsr = compute_lfsr(NegativeProb,ZeroProb)  
    PositiveProb = 1 - NegativeProb - ZeroProb
    lfdr = ZeroProb
    qvalue = qval.from.lfdr(lfdr)
  }
  
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
    if(outputlevel>0){PosteriorMean = PosteriorMean + nonzeromode.fit$nonzeromode}     
  }	   
  
  if(model=="ET"){
    betahat=betahat*sebetahat.orig
    sebetahat = sebetahat.orig
    PosteriorMean = PosteriorMean * sebetahat
    PosteriorSD= PosteriorSD * sebetahat
  }
  
  loglik = calc_gloglik(pi.fit$g, betahat, sebetahat,df, model) 
  
  ##5. Returning the result
  
  result = list(fitted.g=pi.fit$g,call=match.call())
  if (outputlevel>0) {result=c(result,list(PosteriorMean = PosteriorMean,PosteriorSD = PosteriorSD,loglik = loglik))}
  if (outputlevel>1) {result=c(result,list(
                PositiveProb = PositiveProb, NegativeProb = NegativeProb, 
                ZeroProb = ZeroProb,lfsr = lfsr,lfdr = lfdr, qvalue = qvalue, 
                 excludeindex = excludeindex,model = model, optmethod =optmethod))}
  if (outputlevel >2) {result = c(result,list(
    data= list(betahat = betahat, sebetahat = sebetahat,df=df),fit=pi.fit))}
  class(result) = "ash"
  return(result)
  
}

initpi = function(k,n,null.comp,randomstart){
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
  return(pi)
}

setprior=function(prior,k,nullweight,null.comp){
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
  return(prior)
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



#The kth element of this vector is the derivative 
#of the loglik for $\pi=(\pi_0,...,1-\pi_0,...)$ with respect to $\pi_0$ at $\pi_0=1$.
gradient = function(matrix_lik){
  n = nrow(matrix_lik)
  grad = n - colSums(matrix_lik/matrix_lik[,1]) 
  return(grad)
}

#' @title estimate mixture proportions of sigmaa by EM algorithm
#'
#' @param betahat (n vector of observations) 
#' @param sebetahat (n vector of standard errors/deviations of observations)
#' @param g the prior distribution for beta (usually estimated from the data
#' @param prior string, or numeric vector indicating Dirichlet prior on mixture proportions (defaults to "uniform", or (1,1...,1); also can be "nullbiased" (nullweight,1,...,1) to put more weight on first component)
#' @param null.comp the position of the null component
#' @param optmethod name of function to use to do optimization
#' @param df appropriate degrees of freedom for (t) distribution of betahat/sebetahat, default is NULL(Gaussian)
#' @param control A list of control parameters for the SQUAREM algorithm, default value is set to be   control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE). 


#' @return A list, including the final loglikelihood, the null loglikelihood, a n by k likelihoodmatrix with (j,k)th element equal to \eqn{f_k(x_j)},and a flag to indicate convergence.
#
#prior gives the parameter of a Dirichlet prior on pi
#(prior is used to encourage results towards smallest value of sigma when
#likelihood is flat)
#VB provides an approach to estimate the approximate posterior distribution
#of mixture proportions of sigmaa by variational Bayes method
#(use Dirichlet prior and approximate Dirichlet posterior)
#if cxx TRUE use cpp version of R function mixEM

estimate_mixprop = function(betahat,sebetahat,g,prior,optmethod=c("mixEM","mixVBEM","cxxMixSquarem","mixIP"),null.comp=1,df=NULL,control=list()){ 
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  optmethod=match.arg(optmethod)
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)
  
  pi_init = g$pi
  if(optmethod=="mixVBEM"){pi_init=NULL}  #for some reason pi_init doesn't work with mixVBEM
  
  k=ncomp(g)
  n = length(betahat)
  controlinput$tol = min(0.1/n,1.e-7) # set convergence criteria to be more stringent for larger samples
  
  if(controlinput$trace==TRUE){tic()}
  
  matrix_lik = t(compdens_conv(g,betahat,sebetahat,df))
  
  # the last of these conditions checks whether the gradient at the null is negative wrt pi0
  # to avoid running the optimization when the global null (pi0=1) is the optimal.
  if(optmethod=="mixVBEM" || max(prior[-1])>1 || min(gradient(matrix_lik)+prior[1]-1,na.rm=TRUE)<0){
    fit=do.call(optmethod,args = list(matrix_lik= matrix_lik, prior=prior, pi_init=pi_init, control=controlinput))
  } else {
    fit = list(converged=TRUE,pihat=c(1,rep(0,k-1)))
  }
  
  if(!fit$converged){
      warning("Optimization failed to converge. Results may be unreliable. Try increasing maxiter and rerunning.")
  }

  pi = fit$pihat     
  penloglik = penloglik(pi,matrix_lik, prior)
  converged = fit$converged
  
  loglik.final =  penloglik(pi,matrix_lik,1) #compute penloglik without penalty
  null.loglik = sum(log(matrix_lik[,null.comp]))  
  g$pi=pi
  if(controlinput$trace==TRUE){toc()}
  
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
  if(sum(!is.na(lfdr))==0){return(rep(NA,length(lfdr)))}
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
