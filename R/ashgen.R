#' @title Adaptive Shrinkage function for general deconvolution problem
#'
#' @description Takes vectors of estimates (betahat) and the 
#'     error distribution parameters and applies shrinkage to them,
#'     using Empirical Bayes methods, to compute shrunk estimates for
#'     beta.
#'
#' @details Suppose betahat=beta+e, where the errors e come from some error
#'     distribution with CDF & PDF provided, 
#'     we use Empirical methods to compute shrunk estimates for beta.
#'
#' @param betahat a p vector of estimates
#' @param likelihood type of error distribution, can be normal ("normal"), 
#'     generalized t ("t"), log-F ("logF") or self-defined distribution
#'     ("self-defined").
#' @param cdfFUN the CDF function of error distribution. 
#'     Defaults for likelihood "normal", "t" and "logF" are "pnorm", 
#'     "ptgen" and "plogf" respectively. 
#' @param pdfFUN the PDF function of error distribution.
#'     Defaults for likelihood "normal", "t" and "logF" are "dnorm", 
#'     "dtgen" and "dlogf" respectively. 
#' @param FUNargs a list of arguments passed to \code{cdfFUN} & \code{pdfFUN}. 
#'      Each argument must be a scalar or a p vector corresponding to \code{betahat}. 
#' @param etruncFUN (optional) the function to compute expectation of the truncated 
#'      error distribution, which will be used for calculating \code{PosteriorMean}. 
#'      The function should take \code{(a,b,...)}
#'      as inputs and return expectation of truncated error distribution on (\code{a,b}) 
#'      with distribution parameters \code{...}. Default method uses numerical integration
#'      to compute the truncated mean. 
#' @inheritParams ash.workhorse
#'
#' @return ashgen returns an object of \code{\link[base]{class}} "ash", a list with some or all of the following elements (determined by outputlevel) \cr
#' \item{fitted.g}{fitted mixture, either a normalmix or unimix}
#' \item{loglik}{log P(D|mle(pi))}
#' \item{logLR}{log[P(D|mle(pi))/P(D|beta==0)]}
#' \item{PosteriorMean}{A vector consisting the posterior mean of beta from the mixture}
#' \item{PosteriorSD}{A vector consisting the corresponding posterior standard deviation}
#' \item{PositiveProb}{A vector of posterior probability that beta is positive}
#' \item{NegativeProb}{A vector of posterior probability that beta is negative}
#' \item{ZeroProb}{A vector of posterior probability that beta is zero}
#' \item{lfsr}{The local false sign rate}
#' \item{lfdr}{A vector of estimated local false discovery rate}
#' \item{qvalue}{A vector of q values}
#' \item{call}{a call in which all of the specified arguments are specified by their full names}
#' \item{excludeindex}{the vector of index of observations with 0 standard error; if none, then returns NULL}
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
#' #This is equivalent to running beta.ash=ash(betahat,sebetahat)
#' beta.ash = ashgen(betahat,likelihood="normal",FUNargs=list(sd=sebetahat))
#' 
#' #Suppose error distribution is Uniform(min=-2,max=2)
#' beta = c(rep(0,100),rnorm(100))
#' betahat = beta+runif(200,min=-2,max=2)
#' beta.ash = ashgen(betahat,"self-defined",cdfFUN="punif",pdfFUN="dunif",FUNargs=list(min=-2,max=2))
#' 
#' #Provide a self-defined function "etruncunif" to compute mean of
#' #uniform distribution Unif(min,max) truncated to range (a,b)
#' etruncunif = function(a,b,min,max){
#'   if(b<=min | a>=max){
#'     return(0)
#'   }else{
#'     return((max(a,min)+min(b,max))/2)
#'   }
#' }
#' beta.ash = ashgen(betahat,"self-defined",cdfFUN="punif",pdfFUN="dunif",
#'                   FUNargs=list(min=-2,max=2),etruncFUN="etruncunif")
ashgen = function(betahat,
                  likelihood = c("normal","t","logF","self-defined"),
                  cdfFUN = NULL,
                  pdfFUN = NULL,
                  FUNargs = NULL,
                  etruncFUN = NULL,
                  method = c("fdr","shrink"),
                  mixcompdist = c("uniform","halfuniform","+uniform","-uniform"),
                  optmethod = c("mixIP","cxxMixSquarem","mixEM","mixVBEM"),
                  randomstart=FALSE,
                  nullweight=10,nonzeromode=FALSE,
                  pointmass = NULL,
                  prior=c("nullbiased","uniform","unit"),
                  mixsd=NULL, gridmult=sqrt(2),
                  outputlevel=2,
                  g=NULL,
                  fixg=FALSE,
                  cxx=NULL,
                  VB=NULL,
                  control=list()
){
  
  ##1.Handling Input Parameters
  
  method      = match.arg(method)
  mixcompdist = match.arg(mixcompdist)
  optmethod   = match.arg(optmethod)
  likelihood  = match.arg(likelihood)
  
  # Capture all arguments into a list
  oldargs = mget(names(formals()), sys.frame(sys.nframe()))
  newargs = process_args_gen(oldargs)
  
  # Assign each argument in returned list to a variable used by the code next
  for (i in 1:length(newargs)) assign(names(newargs)[i], newargs[[i]])

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
      mixsd = autoselect.mixsd_gen(betahat[completeobs],cdfFUN,pdfFUN, sublist(FUNargs,completeobs), gridmult)
    }
    if(pointmass){
      mixsd = c(0,mixsd)
    }
    
    
    null.comp = which.min(mixsd) #which component is the "null"
    
    k = length(mixsd)
    prior = setprior(prior,k,nullweight,null.comp)
    pi = initpi(k,n,null.comp,randomstart)
    
    if(!is.element(mixcompdist,c("uniform","halfuniform","+uniform","-uniform"))) stop("Error: invalid type of mixcompdist")
    if(mixcompdist=="uniform") g=unimix(pi,-mixsd,mixsd)
    if(mixcompdist=="+uniform") g = unimix(pi,rep(0,k),mixsd)
    if(mixcompdist=="-uniform") g = unimix(pi,-mixsd,rep(0,k))
    if(mixcompdist=="halfuniform"){
      if(min(mixsd)>0){ #simply reflect the components
        g = unimix(c(pi,pi)/2,c(-mixsd,rep(0,k)),c(rep(0,k),mixsd))
        prior = rep(prior, 2)
        pi = rep(pi, 2)
      } else { #define two sets of components, but don't duplicate null component
        null.comp = which.min(mixsd)
        g = unimix(c(pi,pi[-null.comp])/2,c(-mixsd,rep(0,k-1)),c(rep(0,k),mixsd[-null.comp]))
        prior = c(prior,prior[-null.comp])
        pi = c(pi,pi[-null.comp])
      }
    }
  }
  
  #check that all prior are >=1 (as otherwise have problems with infinite penalty)
  if(!all(prior>=1) & optmethod != "mixVBEM"){
    stop("Error: prior must all be >=1 (unless using optmethod mixVBEM)")}
  
  ##3. Fitting the mixture
  if(!fixg){
    pi.fit=estimate_mixprop_gen(betahat[completeobs],cdfFUN,pdfFUN,sublist(FUNargs,completeobs),g,prior,null.comp=null.comp,
                                optmethod=optmethod,control=controlinput)
  } else {
    pi.fit = list(g=g)
  }
  
  ##4. Computing the posterior
  
  n = length(betahat)
  
  PosteriorMean = rep(0,length = n)
  PosteriorSD = rep(0,length = n) 
  
  PosteriorMean[completeobs] = postmean_gen(pi.fit$g,betahat[completeobs],cdfFUN,pdfFUN,sublist(FUNargs,completeobs),
                                            etruncFUN)
  if (sum(is.na(PosteriorMean))){
    warning("NAs produced by posterior mean calculation. Providing a refined etruncFUN may help.")
  }
  #PosteriorSD[completeobs] = postsd(pi.fit$g,betahat[completeobs],cdfFUN,pdfFUN,sublist(FUNargs,completeobs))
  PosteriorSD = rep(NA, length(betahat)) # not implemented for unimix prior
  #FOR MISSING OBSERVATIONS, USE THE PRIOR INSTEAD OF THE POSTERIOR
  PosteriorMean[!completeobs] = calc_mixmean(pi.fit$g)
  #PosteriorSD[!completeobs] = calc_mixsd(pi.fit$g)
  
  ZeroProb = rep(0,length = n)
  NegativeProb = rep(0,length = n)
  ZeroProb[completeobs] = colSums(comppostprob_gen(pi.fit$g,betahat[completeobs],cdfFUN,pdfFUN,sublist(FUNargs,completeobs))
                                  [comp_sd(pi.fit$g) ==0,,drop = FALSE])
  NegativeProb[completeobs] = cdf_post_gen(pi.fit$g, 0, betahat[completeobs],cdfFUN,pdfFUN,sublist(FUNargs,completeobs)) - ZeroProb[completeobs]
  #FOR MISSING OBSERVATIONS, USE THE PRIOR INSTEAD OF THE POSTERIOR
  ZeroProb[!completeobs] = sum(mixprop(pi.fit$g)[comp_sd(pi.fit$g) == 0])
  NegativeProb[!completeobs] = mixcdf(pi.fit$g,0)
  lfsr = compute_lfsr(NegativeProb,ZeroProb)
  PositiveProb = 1 - NegativeProb - ZeroProb
  PositiveProb = ifelse(PositiveProb<0,0,PositiveProb) #deal with numerical issues that lead to numbers <0
  lfdr = ZeroProb
  qvalue = qval.from.lfdr(lfdr)
  svalue = qval.from.lfdr(lfsr)
  
  # kk = ncomp(pi.fit$g)
  # comp_postprob = matrix(0,nrow = kk, ncol = n)
  # comp_postmean = matrix(0,nrow = kk, ncol = n)
  # comp_postmean2 =  matrix(0,nrow = kk, ncol = n)
  # 
  # comp_postprob[,completeobs] = comppostprob(pi.fit$g,betahat[completeobs],sebetahat[completeobs],df)
  # comp_postmean[,completeobs] = comp_postmean(pi.fit$g,betahat[completeobs],sebetahat[completeobs],df)
  # comp_postmean2[,completeobs] = comp_postmean2(pi.fit$g,betahat[completeobs],sebetahat[completeobs],df)
  # 
  # comp_postprob[,!completeobs] = mixprop(pi.fit$g)
  # comp_postmean[,!completeobs] = comp_mean(pi.fit$g)
  # comp_postmean2[,!completeobs] = comp_mean2(pi.fit$g)
  
  loglik = loglik_conv_gen(pi.fit$g,betahat[completeobs],cdfFUN,pdfFUN,sublist(FUNargs,completeobs))
  logLR = loglik-loglik_conv_gen(unimix(1,0,0),betahat[completeobs],cdfFUN,pdfFUN,sublist(FUNargs,completeobs))
  
  ##5. Returning the result
  
  result = list(fitted.g=pi.fit$g,call=match.call(),
                PosteriorMean = PosteriorMean,PosteriorSD = PosteriorSD,
                loglik = loglik, logLR=logLR,
                PositiveProb = PositiveProb, NegativeProb = NegativeProb,
                ZeroProb = ZeroProb,lfsr = lfsr,lfdr = lfdr, qvalue = qvalue, svalue=svalue,
                excludeindex = excludeindex,
                optmethod =optmethod,
                data= list(betahat = betahat, cdfFUN=cdfFUN, pdfFUN=pdfFUN, FUNargs=FUNargs))
  class(result) = "ash"
  return(result)
  
}

# avoid "no visible binding for global variable" note in CRAN check
# These variables are actually defined in process_args
if(getRversion() >= "2.15.1") utils::globalVariables(c("VB","cxx","method","model","mixcompdist","gridmult","control"))

#' Process input arguments for ashgen
#'
#' @param oldargs captured argument list
#' @return list containing the processed arguments
#' @importFrom utils modifyList
process_args_gen = function (oldargs) {
  # Assign each captured argument in the list to a variable
  for (i in 1L:length(oldargs)) assign(names(oldargs)[i], oldargs[[i]])
  
  # Start processing arguments
  # set pdf/cdf functions used for likelihood
  if(likelihood=="normal"){
    cdfFUN = "pnorm"
    pdfFUN = "dnorm"
    if (is.null(FUNargs$mean)){
      FUNargs$mean = 0
    }
    if (is.null(FUNargs$sd)){
      FUNargs$sd = 1
    }
  }else if(likelihood=="t"){
    # use scaled t-distribution
    cdfFUN = "ptgen"
    pdfFUN = "dtgen"
    if (is.null(FUNargs$mean)){
      FUNargs$mean = 0
    }
    if (is.null(FUNargs$sd)){
      FUNargs$sd = 1
    }
    if (is.null(FUNargs$df)){
      stop("df is required for t likelihood!")
    }
  }else if(likelihood=="logF"){
    cdfFUN = "plogf"
    pdfFUN = "dlogf"
    if (is.null(FUNargs$df1) | is.null(FUNargs$df2)){
      stop("df1 & df2 are required for t likelihood!")
    }
  }else if(likelihood=="self-defined"){
    if(!exists(cdfFUN) | !exists(pdfFUN)){
      stop("Error: the self-defined cdf/pdf function does not exist!")
    }
  }
  
  if(!is.null(etruncFUN) & likelihood %in% c("normal","t","logF")){
    warning("Specification of etruncFUN overrides defaults for likelihood normal, t or logF")
  }
  
  # check dimensions of FUNargs
  FUNargs = lapply(FUNargs,dimcheck,n=length(betahat))
  
  # Set optimization method (optmethod)
  
  # if user tries to set both optmethod and VB/cxx that's an error
  if ( !is.null(optmethod) && (!is.null(VB) || !is.null(cxx)) )
    stop("VB and cxx options are deprecated and incompatible with optmethod; use optmethod instead")
  
  # Fallbacks for optmethod
  # By default it will be "mixIP", if REBayes not present then fallback to EM
  if (!requireNamespace("REBayes", quietly = TRUE)) {  # check whether REBayes package is present
    # If REBayes package missing
    message("Due to absence of package REBayes, switching to EM algorithm")
    if (requireNamespace("Rcpp")) {
      optmethod = "cxxMixSquarem"
    } else {
      optmethod = "mixEM"  # fallback if neither Rcpp or REBayes are installed
      message("Using vanilla EM; for faster performance install REBayes (preferred) or Rcpp")
    }
  }
  
  # Check if VB and cxx are set to logical; for backwards compatibility
  if (!is.null(VB)) {
    warning("VB option is deprecated, use optmethod instead")
    if (VB == TRUE) optmethod = "mixVBEM"
  }
  
  if (!is.null(cxx)) {
    warning("cxx option is deprecated, use optmethod instead")
    if (cxx == TRUE)  optmethod = "cxxMixSquarem"
    if (cxx == FALSE) optmethod = "mixEM"
  }
  
  if (optmethod == "mixIP") assertthat::assert_that(requireNamespace("REBayes", quietly = TRUE))
  if (optmethod == "cxxMixSquarem") assertthat::assert_that(requireNamespace("Rcpp", quietly = TRUE))
  
  # method provides a convenient interface to set a particular combinations of parameters for prior an
  # If method is supplied, use it to set up specific values for these parameters; provide warning if values
  # are also specified by user
  # If method is supplied use the user-supplied values (or defaults if user does not specify them)
  if (method == "shrink") {
    
    # Almost equivalent to is.missing(prior)
    if (identical(sort(prior), sort(c("nullbiased","uniform","unit")))) {
      prior = "uniform"
    } else {
      warning("Specification of prior overrides default for method shrink")
    }
    
    if (is.null(pointmass)) {
      pointmass = FALSE
    } else if (pointmass != FALSE) {
      warning("Specification of pointmass overrides default for method shrink")
    }
    
  }
  
  if (method == "fdr") {
    # Almost equivalent to is.missing(prior)
    if (identical(sort(prior), sort(c("nullbiased","uniform","unit")))) {
      prior = "nullbiased"
    } else {
      warning("Specification of prior overrides default for method fdr")
    }
    
    if (is.null(pointmass)) {
      pointmass = TRUE
    } else if (pointmass != TRUE) {
      warning("Specification of pointmass overrides default for method fdr")
    }
    
  }
  
  # Dealing with precise input of betahat, currently we exclude them from the EM algorithm
  betahat.input = betahat
  FUNargs.input = FUNargs
  if (likelihood %in% c("normal","t")){
    excludeindex = c(1:length(FUNargs.input$sd))[FUNargs.input$sd==0]
    if(length(excludeindex) == 0) excludeindex = NULL
    betahat = betahat.input[FUNargs.input$sd != 0]
    FUNargs$sd = FUNargs.input$sd[FUNargs.input$sd != 0]
  }else{
    excludeindex = NULL
  }
  
  # Set observations with infinite standard errors to missing
  # later these missing observations will be ignored in EM, and posterior will be same as prior.
  #sebetahat[sebetahat == Inf] = NA
  #betahat[sebetahat == Inf] = NA
  
  if (mixcompdist == "normal" & !is.null(df))
    stop("Error: Normal mixture for student-t likelihood is not yet implemented")
  
  if (identical(prior, "unit") & optmethod != "mixVBEM")
    stop("Error: unit prior only valid for mixVBEM")
  
  if (mixcompdist == "halfuniform" & !identical(prior, "nullbiased"))
    warning("Use of halfuniform without nullbiased prior can lead to misleading local false sign rates, and so is not recommended")
  
  if (gridmult <= 1) stop("gridmult must be > 1")
  
  #completeobs = (!is.na(betahat) & !is.na(sebetahat))
  completeobs = (!is.na(betahat) & 
                   !is.na(colSums(do.call(rbind,lapply(FUNargs,is.na)))))
  n = sum(completeobs)
  
  # Handling control variables
  control.default = list(K = 1, method = 3, square = TRUE,
                         step.min0 = 1, step.max0 = 1, mstep = 4,
                         kr = 1, objfn.inc = 1, tol = 1.e-07, maxiter = 5000,
                         trace = FALSE)
  if (n > 50000) control.default$trace = TRUE
  namc = names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput = modifyList(control.default, control)
  if (controlinput$maxiter == 0)
    stop("option control$maxiter=0 deprecated; used fixg=TRUE instead")
  
  if (n == 0) stop("Error: all input values are missing")
  
  # Collect everything into a new list
  newargs_names = setdiff(ls(), c("oldargs", "i", "call"))
  newargs = list()
  # assigning NULL to list component will remove that component; so use lapply
  safe_assign = function(x) if (is.null(get(x))) return(NULL) else return(get(x))
  newargs = lapply(newargs_names, safe_assign)
  names(newargs) = newargs_names
  
  return(newargs)
  
}

# check dimensions of list x
# each element of list must be an n-vector or a scalar
dimcheck = function(x,n){
  if (length(x) == 1L){
    x = rep(x, n)
  }
  if (length(x) != n){
    stop("Error: elements in FUNargs must have length 1, or same length as betahat")
  }
  # set Inf to NA
  x[x==Inf] = NA
  return(x)
}

# select sublist from a list
# for each element in this list, select the sub-vector with index idx
sublist = function(list,idx){
  lapply(list, function(x,idx){x[idx]}, idx)
}

#' Estimate mixture proportions of sigmaa by EM algorithm
#' 
#' @inheritParams ashgen
#' @param null.comp the position of the null component
#' @return A list, including the final loglikelihood, the null loglikelihood, 
#'        a n by k likelihoodmatrix with (j,k)th element equal to \eqn{f_k(x_j)},
#'        and a flag to indicate convergence.
estimate_mixprop_gen = function(betahat,cdfFUN,pdfFUN,FUNargs,g,prior,optmethod=c("mixEM","mixVBEM","cxxMixSquarem","mixIP"),null.comp=1,control=list()){
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
  
  matrix_llik = t(log_compdens_conv_gen(g,betahat,cdfFUN,pdfFUN,FUNargs)) #an n by k matrix
  matrix_llik[is.na(matrix_llik)] = log(t(compdens_conv_gen(g,betahat,cdfFUN,pdfFUN,FUNargs)))[is.na(matrix_llik)]
  matrix_llik = matrix_llik - apply(matrix_llik,1, max) #avoid numerical issues by subtracting max of each row
  matrix_lik = exp(matrix_llik)
  
  
  # the last of these conditions checks whether the gradient at the null is negative wrt pi0
  # to avoid running the optimization when the global null (pi0=1) is the optimal.
  if(optmethod=="mixVBEM" || max(prior[-1])>1 || min(gradient_gen(matrix_lik)+prior[1]-1,na.rm=TRUE)<0){
    fit=do.call(optmethod,args = list(matrix_lik= matrix_lik, prior=prior, pi_init=pi_init, control=controlinput))
  } else {
    fit = list(converged=TRUE,pihat=c(1,rep(0,k-1)))
  }
  
  ## check if IP method returns negative mixing proportions. If so, run EM.
  if (optmethod == "mixIP" & (min(fit$pihat) < -10 ^ -12)) {
    message("Interior point method returned negative mixing proportions.\n Switching to EM optimization.")
    optmethod <- "mixEM"
    fit = do.call(optmethod, args = list(matrix_lik = matrix_lik,
                                         prior = prior, pi_init = pi_init,
                                         control = controlinput))
  }
  
  if(!fit$converged){
    warning("Optimization failed to converge. Results may be unreliable. Try increasing maxiter and rerunning.")
  }
  
  pi = fit$pihat
  converged = fit$converged
  niter = fit$niter
  
  loglik.final =  penloglik(pi,matrix_lik,1) #compute penloglik without penalty
  null.loglik = sum(log(matrix_lik[,null.comp]))
  g$pi=pi
  if(controlinput$trace==TRUE){toc()}
  
  return(list(loglik=loglik.final,null.loglik=null.loglik,
              matrix_lik=matrix_lik,converged=converged,g=g,niter=niter))
}

#' @title compdens_conv_gen
#' @description compute the density of the components of the mixture m
#'     when convoluted with a distribution defined by cdfFUN & pdfFUN
#'     with parameters FUNargs, the density is evaluated at x
#' @param m mixture distribution
#' @param x an n vector
#' @inheritParams ashgen
#' @return a k (# of mixture components) by n (length of x) matrix of densities
compdens_conv_gen = function(m,x,cdfFUN,pdfFUN,FUNargs){
  b = pmax(m$b,m$a) #ensure a<b
  a = pmin(m$b,m$a)
  compdens = t(do.call(cdfFUN,c(list(outer(x,a,FUN="-")),FUNargs))-
                 do.call(cdfFUN,c(list(outer(x,b,FUN="-")),FUNargs)))/(b-a)
  compdens[a==b,] = t(do.call(pdfFUN,c(list(outer(x,a,FUN="-")),FUNargs)))[a==b,]
  return(compdens)
}

#' @title log_compdens_conv_gen
#' @description compute the log density of the components of the mixture m
#'     when convoluted with a distribution defined by cdfFUN & pdfFUN
#'     with parameters FUNargs, the density is evaluated at x
#' @inheritParams compdens_conv_gen
#' @return a k (# of mixture components) by n (length of x) matrix of log densities
log_compdens_conv_gen = function(m,x,cdfFUN,pdfFUN,FUNargs){
  b = pmax(m$b,m$a) #ensure a<b
  a = pmin(m$b,m$a)
  lpa = log(do.call(cdfFUN, c(list(outer(x,a,FUN="-")), FUNargs)))
  lpb = log(do.call(cdfFUN, c(list(outer(x,b,FUN="-")), FUNargs)))
  
  lcompdens = t(lpa + log(1-exp(lpb-lpa))) - log(b-a)
  lcompdens[a==b,] = t(log(do.call(pdfFUN, c(list(outer(x,b,FUN="-")),FUNargs))))[a==b,]
  return(lcompdens)
}

#' @title dens_conv_gen
#' @description compute density of mixture m convoluted with a distribution 
#'         defined by cdfFUN & pdfFUN with parameters FUNargs at locations x
#' @inheritParams compdens_conv_gen
dens_conv_gen = function(m,x,cdfFUN,pdfFUN,FUNargs){
  colSums(m$pi * exp(log_compdens_conv_gen(m,x,cdfFUN,pdfFUN,FUNargs)))
}


#' @title loglik_conv_gen
#' @description find log likelihood of data in betahat, when the
#'     mixture m is convolved with a distribution defined
#'     by cdfFUN & pdfFUN with parameters FUNargs
#' @inheritParams compdens_conv_gen
loglik_conv_gen = function(m,x,cdfFUN,pdfFUN,FUNargs){
  sum(log(dens_conv_gen(m,x,cdfFUN,pdfFUN,FUNargs)))
}

#The kth element of this vector is the derivative
#of the loglik for $\pi=(\pi_0,...,1-\pi_0,...)$ with 
#respect to $\pi_0$ at $\pi_0=1$.
gradient_gen = function(matrix_lik){
  n = nrow(matrix_lik)
  # avoid dividing by 0
  grad = n - colSums(matrix_lik/(matrix_lik[,1]+1e-200))
  return(grad)
}

#' @title postmean_gen
#' @description output posterior mean for beta for prior mixture m,
#' given observations betahat
#' @param m mixture distribution with k components
#' @inheritParams ashgen
postmean_gen = function(m,betahat,cdfFUN,pdfFUN,FUNargs,etruncFUN){
  colSums(comppostprob_gen(m,betahat,cdfFUN,pdfFUN,FUNargs) * 
            comp_postmean_gen(m,betahat,cdfFUN,pdfFUN,FUNargs,etruncFUN))
}

#' @title comppostprob_gen
#' @description compute the posterior prob that each observation came
#'     from each component of the mixture m, output a k by n vector of
#'     probabilities computed by weighting the component densities by
#'     pi and then normalizing
#' @inheritParams compdens_conv_gen
comppostprob_gen = function(m,x,cdfFUN,pdfFUN,FUNargs){
  lpost = log_compdens_conv_gen(m,x,cdfFUN,pdfFUN,FUNargs) + log(m$pi) # lpost is k by n of log posterior prob (unnormalized)
  lpost[is.na(lpost)] = (log(compdens_conv_gen(m,x,cdfFUN,pdfFUN,FUNargs)) + log(m$pi))[is.na(lpost)]
  lpmax = apply(lpost,2,max) #dmax is of length n
  tmp = exp(t(lpost)-lpmax) #subtracting the max of the logs is just done for numerical stability
  tmp = tmp/rowSums(tmp)
  #ismissing = (is.na(x) | is.na(s))
  ismissing = (is.na(x) | is.na(colSums(do.call(rbind,lapply(FUNargs,is.na)))))
  tmp[ismissing,]=m$pi
  t(tmp)
}

#' @title comp_postmean
#' @description output posterior mean for beta for each component of
#'     prior mixture m, given observations betahat and FUNargs
#' @inheritParams postmean_gen
comp_postmean_gen = function(m,betahat,cdfFUN,pdfFUN,FUNargs,etruncFUN){
  # if type==normal,t,logf, just use the refined my_etruncnorm....
  if((cdfFUN=="pnorm" & pdfFUN=="dnorm") |
     (cdfFUN=="ptgen" & pdfFUN=="dtgen")){
    comp_postmean(m,betahat,FUNargs$sd,FUNargs$df)
  }else if(cdfFUN=="plogf" & pdfFUN=="dlogf"){
    #comp_postmean_logf(m,betahat,FUNargs$df1[1],FUNargs$df2[1])
    comp_postmean_logf(m,betahat,FUNargs$df1,FUNargs$df2)
  }else{
    alpha = outer(betahat, m$b,FUN="-")
    beta = outer(betahat, m$a, FUN="-")
    tmp = betahat - do.call(my_etruncdistn, 
                           list(a=alpha, b=beta, cdfFUN=cdfFUN, pdfFUN=pdfFUN,
                                FUNargs=FUNargs, etruncFUN=etruncFUN))
    ismissing = (is.na(betahat) | is.na(colSums(do.call(rbind,lapply(FUNargs,is.na)))))
    tmp[ismissing,]= (m$a+m$b)/2
    t(tmp)
  }
}

#' @title my_etruncdistn
#' @description Compute expectation of truncated error distribution.
#'
#' @param a Left limit of distribution.
#' @param b Right limit of distribution.
#' @inheritParams ashgen
my_etruncdistn = function(a,b,cdfFUN,pdfFUN,FUNargs,etruncFUN){
  if(!is.null(etruncFUN)){
    tmp = do.call(mapply, c(list(FUN=etruncFUN,a=c(a),b=c(b)), FUNargs))
  }else{
    tmp = do.call(mapply, c(list(FUN=my_etruncdistn_single,a=c(a),b=c(b),
                                 cdfFUN=cdfFUN,pdfFUN=pdfFUN), FUNargs))
  }
  tmp = matrix(tmp, nrow=dim(a)[1])
  return(tmp)
}

# compute expectation of truncated error distribution
# for scalars a and b
my_etruncdistn_single = function(a,b,cdfFUN,pdfFUN,...){
  if(a == b){
    tmp = a
  }else{
    # numerical integration
    xpdf = function(x,pdfFUN,...){
      x*do.call(pdfFUN,list(x,...))
    }
    tmp = try(integrate(xpdf,a,b,pdfFUN,...)$value,silent=TRUE)
    
    if (class(tmp)!="try-error"){
      # no probability mass on (a,b), expectation is 0
      if (do.call(cdfFUN,list(b,...))
          -do.call(cdfFUN,list(a,...))==0){
        tmp = 0
      }else{
        tmp = tmp/(do.call(cdfFUN,list(b,...))
                   -do.call(cdfFUN,list(a,...)))
      }
    }else{
      tmp = NA
    }
  }
  return(tmp)
}

#' @title compcdf_post_gen
#' @description evaluate cdf of posterior distribution of beta at c. m
#'     is the prior on beta, a mixture; c is location of evaluation
#'     assumption is betahat | beta = beta+e, where e come from error 
#'     distribution defined by cdfFUN, pdfFUN and FUNargs.
#' @param m mixture distribution with k components
#' @param c a scalar
#' @inheritParams ashgen
compcdf_post_gen=function(m,c,betahat,cdfFUN,pdfFUN,FUNargs){
  k = length(m$pi)
  n=length(betahat)
  tmp = matrix(1,nrow=k,ncol=n)
  tmp[m$a > c,] = 0
  subset = m$a<=c & m$b>c # subset of components (1..k) with nontrivial cdf
  if(sum(subset)>0){
    pna = do.call(cdfFUN, c(list(outer(betahat,m$a[subset],FUN="-")), FUNargs))
    pnc = do.call(cdfFUN, c(list(outer(betahat,rep(c,sum(subset)),FUN="-")), FUNargs))
    pnb = do.call(cdfFUN, c(list(outer(betahat,m$b[subset],FUN="-")), FUNargs))
    
    tmp[subset,] = t((pnc-pna)/(pnb-pna))
  }
  subset = (m$a == m$b) #subset of components with trivial cdf
  tmp[subset,]= rep(m$a[subset] <= c,n)
  #Occasionally we would encounter issue such that in some entries pna[i,j]=pnb[i,j]=pnc[i,j]=0 or pna=pnb=pnc=1
  #Those are the observations with significant betahat(small sebetahat), resulting in pnorm() return 1 or 0
  #due to the thin tail property of normal distribution.(or t-distribution, although less likely to occur)
  #Then R would be dividing 0 by 0, resulting  in NA values
  #In practice, those observations would have 0 probability of belonging to those "problematic" components
  #Thus any sensible value in [0,1] would not matter much, as they are highly unlikely to come from those
  #components in posterior distribution.
  #Here we simply assign the "naive" value as as (c-a)/(b-a)
  #As the component pdf is rather smaller over the region.
  tmpnaive=matrix(rep((c-m$a)/(m$b-m$a),length(betahat)),nrow=k,ncol=n)
  tmp[is.nan(tmp)]= tmpnaive[is.nan(tmp)]
  tmp
}

#' @title cdf_post_gen
#' @description evaluate cdf of posterior distribution of beta at c. m
#'     is the prior on beta, a mixture; c is location of evaluation
#'     assumption is betahat | beta = beta+e, where e come from error 
#'     distribution defined by cdfFUN, pdfFUN and FUNargs.
#' @inheritParams compcdf_post_gen
cdf_post_gen=function(m,c,betahat,cdfFUN,pdfFUN,FUNargs){
  colSums(comppostprob_gen(m,betahat,cdfFUN,pdfFUN,FUNargs)*
            compcdf_post_gen(m,c,betahat,cdfFUN,pdfFUN,FUNargs))
}

# try to select a default range for the sigmaa values
# that should be used, based on the values of betahat and sebetahat
# mult is the multiplier by which the sds differ across the grid
autoselect.mixsd_gen = function(betahat,cdfFUN,pdfFUN,FUNargs,mult){
  if((cdfFUN=="pnorm" & pdfFUN=="dnorm") |
     (cdfFUN=="ptgen" & pdfFUN=="dtgen")){
    return(autoselect.mixsd(betahat, FUNargs$sd, mult))
  }else if(cdfFUN=="plogf" & pdfFUN=="dlogf"){
    return(autoselect.mixsd_logf(betahat, FUNargs$df1, FUNargs$df2, mult))
  }else{
    sd = rep(1, length(betahat))
    for (i in 1:length(betahat)){
      currargs = lapply(FUNargs,`[[`,i)
      sd[i] = distn_sd(pdfFUN, currargs, range=min(2*max(abs(betahat)),10))
    }
    sd[sd==0] = sd(betahat)
    return(autoselect.mixsd(betahat, sd, mult))
  }
}

# compute standard deviation of a distribution
# pdfFUN: pdf function
# FUNargs: arguments for pdfFUN, each element is scalar
# range: do numerical integration on (-range, range)
distn_sd = function(pdfFUN, FUNargs, range=10){
  xpdf = function(x,pdfFUN,FUNargs){
    x*do.call(pdfFUN,c(list(x),FUNargs))
  }
  x2pdf = function(x,pdfFUN,FUNargs){
    x^2*do.call(pdfFUN,c(list(x),FUNargs))
  }
  # second moment
  moment2 = try(integrate(x2pdf,-range,range,pdfFUN,FUNargs)$value, silent=TRUE)
  # first moment
  moment1 = try(integrate(xpdf,-range,range,pdfFUN,FUNargs)$value, silent=TRUE)
  if (class(moment1)=="numeric" & class(moment2)=="numeric"){
    return(sqrt(moment2-moment1^2))
  }else{
    return(NA)
  }
}

#' @title Generalized student t distribution
#' @description Distribution function for the generalized t distribution with 
#'       \code{df} degrees of freedom, location parameter \code{mean}, scale parameter \code{sd} 
#'       (and optional non-centrality parameter \code{ncp}). 
#'       i.e. (X-\code{mean})/\code{sd} follows standard t-distribution with \code{df} and \code{ncp}.
#' @param q vector of quantiles
#' @param df degrees of freedom (> 0, maybe non-integer). 
#' @param ncp non-centrality parameter delta; only for abs(ncp) <= 37.62. 
#'        If omitted, use the central t distribution.
#' @param mean location parameter
#' @param sd scale parameter
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X ≤ x], otherwise, P[X > x].
#' @return The distribution function.
#' @export
ptgen = function(q, df, ncp, mean=0, sd=1, lower.tail=TRUE, log.p=FALSE){
  return(pt((q-mean)/sd, df=df, ncp=ncp, lower.tail=lower.tail, log.p=log.p))
}

#' @title Generalized student t distribution
#' @description Density function for the generalized t distribution with 
#'       \code{df} degrees of freedom, location parameter \code{mean}, scale parameter \code{sd} 
#'       (and optional non-centrality parameter \code{ncp}). 
#'       i.e. (X-\code{mean})/\code{sd} follows standard t-distribution with \code{df} and \code{ncp}.
#' @param x vector of quantiles
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @inheritParams ptgen
#' @return The density.
#' @export
dtgen = function(x, df, ncp, mean=0, sd=1, log=FALSE){
  return(dt((x-mean)/sd, df=df, ncp=ncp, log=log))
}

#' @title The log-F distribution
#' @description Distribution function for the log-F distribution with \code{df1} and \code{df2}
#'  degrees of freedom (and optional non-centrality parameter \code{ncp}).
#' @param q vector of quantiles
#' @param df1,df2 degrees of freedom
#' @param ncp non-centrality parameter. If omitted the central F is assumed.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X ≤ x], otherwise, P[X > x].
#' @return The distribution function. 
#' @export
plogf = function(q, df1, df2, ncp, lower.tail=TRUE, log.p=FALSE){
  return(pf(exp(q), df1=df1, df2=df2, ncp=ncp, lower.tail=lower.tail,
            log.p=log.p))
}

#' @title The log-F distribution
#' @description Density function for the log-F distribution with \code{df1} and \code{df2}
#'  degrees of freedom (and optional non-centrality parameter \code{ncp}).
#' @param x vector of quantiles
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @inheritParams plogf
#' @return The density function. 
#' @export
dlogf = function(x, df1, df2, ncp, log=FALSE){
  if (log==FALSE){
    df(exp(x), df1=df1, df2=df2, ncp=ncp)*exp(x)
  }else{
    df(exp(x), df1=df1, df2=df2, ncp=ncp, log=TRUE)+x
  }
}

#posterior mean for alpha for prior mixture m,
#given observations logfhat=alpha+e where 
#e~logF(df1,df2)
postmean_logf = function(m,logfhat,v1,v2){
  colSums(comppostprob_logf(m,logfhat,v1,v2) * comp_postmean_logf(m,logfhat,v1,v2))
}

#posterior component-wise mean for alpha for prior mixture m,
#given observations logfhat=alpha+e where 
#e~logF(df1,df2)
comp_postmean_logf = function(m,logfhat,v1,v2){
  alpha = outer(-logfhat, m$a,FUN="+")
  beta = outer(-logfhat, m$b, FUN="+")
  tmp = matrix(mapply(my_etrunclogf,c(alpha),c(beta),v2,v1),
               nrow=length(logfhat))
  tmp = logfhat+tmp
  
  ismissing = is.na(logfhat)
  tmp[ismissing,]= (m$a+m$b)/2
  t(tmp)
}

#' @title my_etruncnorm
#' @description Compute expectation of truncated log-F distribution.
#'
#' @param a Left limit of distribution.
#' @param b Right limit of distribution.
#' @param df1,df2 degrees of freedom
#' @export
my_etrunclogf= function(a,b,df1,df2){
  if (a==b){
    tmp = a
  }else{
    tmp = try(etrunclogf(df1=df1, df2=df2, a=a, b=b, adj=FALSE),silent=TRUE)
    if (class(tmp)=="try-error"){
      tmp = try(etrunclogf(df1=df1, df2=df2, a=a, b=b, adj=TRUE),silent=TRUE)
    }
    
    if (class(tmp)=="try-error"){
      #tmp = NA
      tmp = (a+b)/2
    }
  }
  return(tmp) #deal with extreme case a=b
}

# x*dlogf
etrunclogf_num = function(x,df1,df2,a,b){
  #multiply c to avoid numerical issues
  c = 10^(-round(min(log10(df(exp(a),df1,df2)*exp(a)),
                     log10(df(exp(b),df1,df2)*exp(b)))))
  c*x*df(exp(x),df1=df1,df2=df2)*exp(x)
}

# dlogf
etrunclogf_denom = function(x,df1,df2,a,b){
  #multiply c to avoid numerical issues
  c = 10^(-round(min(log10(df(exp(a),df1,df2)*exp(a)),
                     log10(df(exp(b),df1,df2)*exp(b)))))
  c*df(exp(x),df1=df1,df2=df2)*exp(x)
}

# x multiply by the density of truncated log-F distribution on (a,b) at x
xdtrunclogf = function(x,df1,df2,a,b){
  x*df(exp(x),df1=df1,df2=df2)*exp(x)/(pf(exp(b),df1,df2)-pf(exp(a),df1,df2))
}

# compute expectation of truncated log-F distribution.
etrunclogf = function(df1,df2,a,b,adj=FALSE){
  if (adj==TRUE){
    # numerator and denominator both multiply a constant to avoid numerical issues
    n = integrate(etrunclogf_num, lower=a, upper=b, df1=df1,df2=df2,a=a,b=b)$value
    d = integrate(etrunclogf_denom, lower=a, upper=b, df1=df1,df2=df2,a=a,b=b)$value
    return(n/d)
  }else{
    return(integrate(xdtrunclogf, lower=a, upper=b, df1=df1,df2=df2,a=a,b=b)$value)
  } 
}

# try to select a default range for the sigmaa values
# that should be used, based on the values of logfhat, df1 and df2
# mult is the multiplier by which the sds differ across the grid
autoselect.mixsd_logf = function(logfhat,df1,df2,mult){
  sigmaamin = log(qf(0.85,df1=df1,df2=df2))/10 #so that the minimum is small compared with measurement precision
  if(all(logfhat^2<=log(qf(0.85,df1,df2))^2)){
    sigmaamax = 8*sigmaamin #to deal with the occassional odd case where this could happen; 8 is arbitrary
  }else{
    sigmaamax = 2*sqrt(max(logfhat^2-log(qf(0.85,df1,df2))^2)) #this computes a rough largest value you'd want to use, based on idea that sigmaamax^2 + sebetahat^2 should be at least betahat^2   
  }
  if(mult==0){
    return(c(0,sigmaamax/2))
  }else{
    npoint = ceiling(log2(sigmaamax/sigmaamin)/log2(mult))
    return(mult^((-npoint):0) * sigmaamax)
  }
}