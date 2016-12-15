#' @useDynLib ashr
#' @import truncnorm SQUAREM doParallel pscl Rcpp foreach parallel

#' @title Main Adaptive Shrinkage function
#'
#' @description Takes vectors of estimates (betahat) and their
#'     standard errors (sebetahat), together with degrees of freedom (df)
#'     and applies shrinkage to them, using Empirical Bayes methods, to compute shrunk estimates for
#'     beta.
#'
#' @details This function is actually just a simple wrapper that
#'     passes its parameters to \code{\link{ash.workhorse}} which
#'     provides more documented options for advanced use. See readme
#'     for more details.
#'
#' @param betahat a p vector of estimates
#' @param sebetahat a p vector of corresponding standard errors
#' @param mixcompdist distribution of components in mixture
#'     ("uniform","halfuniform" or "normal"; "+uniform" or
#'     "-uniform"), the default is "uniform". If you believe your
#'     effects may be asymmetric, use "halfuniform". If you want to
#'     allow only positive/negative effects use "+uniform"/"-uniform".
#'     The use of "normal" is permitted only if df=NULL.
#' @param df appropriate degrees of freedom for (t) distribution of
#'     betahat/sebetahat, default is NULL which is actually treated as
#'     infinity (Gaussian)
#' @param ... Further arguments to be passed to
#'     \code{\link{ash.workhorse}}.
#'
#' @return ash returns an object of \code{\link[base]{class}} "ash", a list with some or all of the following elements (determined by outputlevel) \cr
#' \item{fitted_g}{fitted mixture}
#' \item{loglik}{log P(D|fitted_g)}
#' \item{logLR}{log[P(D|fitted_g)/P(D|beta==0)]}
#' \item{result}{A dataframe whose columns are}
#' \describe{
#'  \item{NegativeProb}{A vector of posterior probability that beta is negative}
#'  \item{PositiveProb}{A vector of posterior probability that beta is positive}
#'  \item{lfsr}{A vector of estimated local false sign rate}
#'  \item{lfdr}{A vector of estimated local false discovery rate}
#'  \item{qvalue}{A vector of q values}
#'  \item{svalue}{A vector of s values}
#'  \item{PosteriorMean}{A vector consisting the posterior mean of beta from the mixture}
#'  \item{PosteriorSD}{A vector consisting the corresponding posterior standard deviation}
#'  }
#' \item{call}{a call in which all of the specified arguments are specified by their full names}
#' \item{data}{a list containing details of the data and models used (mostly for internal use)}
#' \item{fit_details}{a list containing results of mixture optimization, and matrix of component log-likelihoods used in this optimization}
#'
#' @seealso \code{\link{ash.workhorse}} for complete specification of ash function
#' @seealso \code{\link{ashci}} for computation of credible intervals after getting the ash object return by \code{ash()}
#'
#' @export
#' @examples
#' beta = c(rep(0,100),rnorm(100))
#' sebetahat = abs(rnorm(200,0,1))
#' betahat = rnorm(200,beta,sebetahat)
#' beta.ash = ash(betahat, sebetahat)
#' names(beta.ash)
#' head(beta.ash$result) # the main dataframe of results
#' graphics::plot(betahat,beta.ash$result$PosteriorMean,xlim=c(-4,4),ylim=c(-4,4))
#'
#' CIMatrix=ashci(beta.ash,level=0.95)
#' print(CIMatrix)
#'
#' #Illustrating the non-zero mode feature
#' betahat=betahat+5
#' beta.ash = ash(betahat, sebetahat)
#' graphics::plot(betahat,beta.ash$result$PosteriorMean)
#' betan.ash=ash(betahat, sebetahat,mode=5)
#' graphics::plot(betahat, betan.ash$result$PosteriorMean)
#' summary(betan.ash)
ash <- function (betahat, sebetahat,
                 mixcompdist = c("uniform","halfuniform","normal","+uniform",
                                 "-uniform"),
                 df = NULL,...)

  # TO DO: Explain here what this does. It certainly isn't clear
  # (thanks in part to R's strangeness)!
  utils::modifyList(ash.workhorse(betahat,sebetahat,
                                  mixcompdist = mixcompdist,df = df,...),
                    list(call = match.call()))

#' @title Detailed Adaptive Shrinkage function
#'
#' @description Takes vectors of estimates (betahat) and their
#'     standard errors (sebetahat), and applies shrinkage to them,
#'     using Empirical Bayes methods, to compute shrunk estimates for
#'     beta. This is the more detailed version of ash for "research"
#'     use.  Most users will be happy with the ash function, which
#'     provides the same usage, but documents only the main options
#'     for simplicity.
#'
#' @details See readme for more details.
#'
#' @param betahat a p vector of estimates
#' @param sebetahat a p vector of corresponding standard errors
#' @param method specifies how ash is to be run. Can be "shrinkage"
#'     (if main aim is shrinkage) or "fdr" (if main aim is to assess
#'     fdr or fsr) This is simply a convenient way to specify certain
#'     combinations of parameters: "shrinkage" sets pointmass=FALSE
#'     and prior="uniform"; "fdr" sets pointmass=TRUE and
#'     prior="nullbiased".
#' @param mixcompdist distribution of components in mixture (
#'     "uniform","halfuniform","normal" or "+uniform"), the default
#'     value is "uniform" use "halfuniform" to allow for assymetric g,
#'     and "+uniform"/"-uniform" to constrain g to be
#'     positive/negative.
#' @param optmethod specifies the function implementing an optimization method. Default is
#'     "mixIP", an interior point method, if REBayes is installed;
#'     otherwise an EM algorithm is used. The interior point method is
#'     faster for large problems (n>2000), particularly when method="shrink".
#' @param df appropriate degrees of freedom for (t) distribution of
#'     betahat/sebetahat, default is NULL(Gaussian)
#' @param nullweight scalar, the weight put on the prior under
#'     "nullbiased" specification, see \code{prior}
#' @param mode either numeric (indicating mode of g) or string "estimate",
#'      to indicate mode should be estimated.
#' @param pointmass logical, indicating whether to use a point mass at
#'     zero as one of components for a mixture distribution
#' @param prior string, or numeric vector indicating Dirichlet prior
#'     on mixture proportions (defaults to "uniform", or (1,1...,1);
#'     also can be "nullbiased" (nullweight,1,...,1) to put more
#'     weight on first component), or "unit" (1/K,...,1/K) [for
#'     optmethod=mixVBEM version only]
#' @param mixsd vector of sds for underlying mixture components
#' @param gridmult the multiplier by which the default grid values for
#'     mixsd differ by one another. (Smaller values produce finer
#'     grids)
#' @param outputlevel determines amount of output. There are several numeric options [0=just fitted g;
#'     1=also PosteriorMean and PosteriorSD; 2= everything usually
#'     needed; 3=also include results of mixture fitting procedure
#'     (includes matrix of log-likelihoods used to fit mixture); 4=
#'     output additional things required by flash (flash_data)]. Otherwise the user can also specify
#'     the output they require in detail (see Examples)
#' @param g the prior distribution for beta (usually estimated from
#'     the data; this is used primarily in simulated data to do
#'     computations with the "true" g)
#' @param fixg if TRUE, don't estimate g but use the specified g -
#'     useful for computations under the "true" g in simulations
#' @param alpha numeric value of alpha parameter in the model
#' @param control A list of control parameters passed to optmethod
#' @param lik contains details of the likelihood used; for general ash
#'
#' @return ash returns an object of \code{\link[base]{class}} "ash", a list with some or all of the following elements (determined by outputlevel) \cr
#' \item{fitted_g}{fitted mixture, either a normalmix or unimix}
#' \item{loglik}{log P(D|mle(pi))}
#' \item{logLR}{log[P(D|mle(pi))/P(D|beta==0)]}
#' \item{result}{A dataframe whose columns are}
#' \describe{
#'  \item{NegativeProb}{A vector of posterior probability that beta is negative}
#'  \item{PositiveProb}{A vector of posterior probability that beta is positive}
#'  \item{lfsr}{A vector of estimated local false sign rate}
#'  \item{lfdr}{A vector of estimated local false discovery rate}
#'  \item{qvalue}{A vector of q values}
#'  \item{svalue}{A vector of s values}
#'  \item{PosteriorMean}{A vector consisting the posterior mean of beta from the mixture}
#'  \item{PosteriorSD}{A vector consisting the corresponding posterior standard deviation}
#'  }
#' \item{call}{a call in which all of the specified arguments are specified by their full names}
#' \item{data}{a list containing details of the data and models used (mostly for internal use)}
#' \item{fit_details}{a list containing results of mixture optimization, and matrix of component log-likelihoods used in this optimization}
#'
#' @seealso \code{\link{ash}} for simplified specification of ash function
#' @seealso \code{\link{ashci}} for computation of credible intervals
#'     after getting the ash object return by \code{ash()}
#'
#' @export
#' @examples
#' beta = c(rep(0,100),rnorm(100))
#' sebetahat = abs(rnorm(200,0,1))
#' betahat = rnorm(200,beta,sebetahat)
#' beta.ash = ash(betahat, sebetahat)
#' names(beta.ash)
#' head(beta.ash$result) #dataframe of results
#' head(get_lfsr(beta.ash)) #get lfsr
#' head(get_pm(beta.ash)) #get posterior mean
#' graphics::plot(betahat,get_pm(beta.ash),xlim=c(-4,4),ylim=c(-4,4))
#'
#' CIMatrix=ashci(beta.ash,level=0.95) #note currently default is only compute CIs for lfsr<0.05
#' print(CIMatrix)
#'
#' #Running ash with different error models
#' beta.ash1 = ash(betahat, sebetahat, lik = normal_lik())
#' beta.ash2 = ash(betahat, sebetahat, lik = t_lik(df=4))
#'
#' e = rnorm(100)+log(rf(100,df1=10,df2=10)) # simulated data with log(F) error
#' e.ash = ash(e,1,lik=logF_lik(df1=10,df2=10))
#'
#' # Specifying the output
#' beta.ash = ash(betahat, sebetahat, output = c("fitted_g","logLR","lfsr"))
#'
#' #Illustrating the non-zero mode feature
#' betahat=betahat+5
#' beta.ash = ash(betahat, sebetahat)
#' graphics::plot(betahat,beta.ash$result$PosteriorMean)
#' betan.ash=ash(betahat, sebetahat,mode=5)
#' graphics::plot(betahat, betan.ash$result$PosteriorMean)
#' summary(betan.ash)
#'
#' #Running ash with a pre-specified g, rather than estimating it
#' beta = c(rep(0,100),rnorm(100))
#' sebetahat = abs(rnorm(200,0,1))
#' betahat = rnorm(200,beta,sebetahat)
#' true_g = normalmix(c(0.5,0.5),c(0,0),c(0,1)) # define true g
#' ## Passing this g into ash causes it to i) take the sd and the means
#' ## for each component from this g, and ii) initialize pi to the value
#' ## from this g.
#' beta.ash = ash(betahat, sebetahat,g=true_g,fixg=TRUE)
ash.workhorse <-
    function(betahat, sebetahat, method = c("fdr","shrink"),
             mixcompdist = c("uniform","halfuniform","normal","+uniform",
                             "-uniform"),
             optmethod = c("mixIP","cxxMixSquarem","mixEM","mixVBEM"),
             df = NULL,nullweight = 10,pointmass = TRUE,
             prior = c("nullbiased","uniform","unit"),mixsd = NULL,
             gridmult = sqrt(2),outputlevel = 2,g = NULL,fixg = FALSE,
             mode = 0,alpha = 0,control = list(),lik = NULL) {

  if(!missing(pointmass) & !missing(method))
    stop("Specify either method or pointmass, not both")
  if(!missing(prior) & !missing(method))
    stop("Specify either method or prior, not both")
  if(!missing(mode) & !missing(g))
    stop("Specify either mode or g, not both")
  if(!missing(method)){
    method = match.arg(method)
    if (method == "shrink"){pointmass =FALSE; prior="uniform"}
    if (method == "fdr"){pointmass =TRUE; prior= "nullbiased"}
  }

  ## Check to see if is Inf, then switch to NULL.
  if (!is.null(df)) {
    if (df == Inf) {
      df <- NULL
    }
  }

  if(mode=="estimate"){ #just pass everything through to ash.estmode for non-zero-mode
    args <- as.list(environment())
    args$mode = NULL
    args$outputlevel = NULL
    args$method=NULL # avoid specifying method as well as prior/pointmass
    args$g = NULL # avoid specifying g as well as mode
    #args = as.list( match.call() )
    mode = do.call(ash.estmode,args)}


  ##1.Handling Input Parameters
  mixcompdist = match.arg(mixcompdist)
  optmethod   = match.arg(optmethod)
  prior       = match.arg(prior)

  # Set optimization method
  optmethod = set_optmethod(optmethod)
  check_args(mixcompdist,df,prior,optmethod,gridmult,sebetahat,betahat)
  if(is.null(lik)){ #set likelihood based on defaults if missing
    if(is.null(df)){
      lik = normal_lik()
    } else {lik = t_lik(df)}
  }
  check_lik(lik) # minimal check that it obeys requirements
  lik = add_etruncFUN(lik) #if missing, add a function to compute mean of truncated distribution
  data = set_data(betahat, sebetahat, lik, alpha)

  ##2. Generating mixture distribution g

  if(fixg & missing(g)){stop("if fixg=TRUE then you must specify g!")}

  if(!is.null(g)){
    k=ncomp(g)
    null.comp=1 #null.comp not actually used
    prior = setprior(prior,k,nullweight,null.comp)
  } else {
    if(is.null(mixsd)){
      mixsd = autoselect.mixsd(data,gridmult,mode)
    }
    if(pointmass){
      mixsd = c(0,mixsd)
    }
    null.comp = which.min(mixsd) #which component is the "null"

    k = length(mixsd)
    prior = setprior(prior,k,nullweight,null.comp)
    pi = initpi(k,length(data$x),null.comp)

    if(!is.element(mixcompdist,c("normal","uniform","halfuniform","+uniform","-uniform")))
      stop("Error: invalid type of mixcompdist")
    if(mixcompdist == "normal") g=normalmix(pi,rep(mode,k),mixsd)
    if(mixcompdist == "uniform") g=unimix(pi,mode - mixsd,mode + mixsd)
    if(mixcompdist == "+uniform") g = unimix(pi,rep(mode,k),mode+mixsd)
    if(mixcompdist == "-uniform") g = unimix(pi,mode-mixsd,rep(mode,k))
    if(mixcompdist == "halfuniform"){
      if(min(mixsd)>0){ #simply reflect the components
        g = unimix(c(pi,pi)/2,c(mode-mixsd,rep(mode,k)),c(rep(mode,k),mode+mixsd))
        prior = rep(prior, 2)
        pi = rep(pi, 2)
      } else { #define two sets of components, but don't duplicate null component
        null.comp=which.min(mixsd)
        g = unimix(c(pi,pi[-null.comp])/2,c(mode-mixsd,rep(mode,k-1)),c(rep(mode,k),mode+mixsd[-null.comp]))
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
    pi.fit=estimate_mixprop(data,g,prior,optmethod=optmethod,control=control)
  } else {
    pi.fit = list(g=g)
  }

  ##4. Computing the return values

  val = list() # val will hold the return value
  ghat = pi.fit$g
  output = set_output(outputlevel) #sets up flags for what to output
  if("fitted_g" %in% output){val = c(val,list(fitted_g=ghat))}
  if("loglik" %in% output){val = c(val,list(loglik =calc_loglik(ghat,data)))}
  if("logLR" %in% output){val = c(val,list(logLR=calc_logLR(ghat,data)))}
  if("data" %in% output){val = c(val,list(data=data))}
  if("fit_details" %in% output){val = c(val,list(fit_details = pi.fit))}
  if("flash_data" %in% output){val = c(val, list(flash_data=calc_flash_data(ghat,data)))}

  # Compute the result component of value -
  # result is a dataframe containing lfsr, etc
  # resfns is a list of functions used to produce columns of that dataframe
  resfns = set_resfns(output)
  if(length(resfns)>0){
    result = data.frame(betahat = betahat,sebetahat = sebetahat)
    if(!is.null(df)){result$df = df}
    result = cbind(result,as.data.frame(lapply(resfns,do.call,list(g=pi.fit$g,data=data))))
    val = c(val, list(result=result))
  }

  ##5. Returning the val

  class(val) = "ash"
  return(val)

}


#adds result of applying f to (g,data) to the list res
#the result of f should be a list with named elements
add_list = function(f,g,data,res){
  return(c(res,do.call(f, list(g=g,data=data))))
}

initpi = function(k,n,null.comp,randomstart=FALSE){
  if(randomstart){
    pi = stats::rgamma(k,1,1)
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
    } else if(prior=="unit"){
      prior = rep(1/k,k)
    }
  }
  if(length(prior)!=k | !is.numeric(prior)){
    stop("invalid prior specification")
  }
  return(prior)
}


#' Function to compute the local false sign rate
#'
#' @param NegativeProb A vector of posterior probability that beta is
#'     negative.
#' @param ZeroProb A vector of posterior probability that beta is
#'     zero.
#' @return The local false sign rate.
compute_lfsr = function(NegativeProb,ZeroProb){
  ifelse(NegativeProb> 0.5*(1-ZeroProb),1-NegativeProb,NegativeProb+ZeroProb)
}



#The kth element of this vector is the derivative
#of the loglik for $\pi=(\pi_0,...,1-\pi_0,...)$ with respect to $\pi_0$ at $\pi_0=1$.
gradient = function(matrix_lik){
  n = nrow(matrix_lik)
  grad = n - colSums(matrix_lik/matrix_lik[,1])
  return(grad)
}

#' Estimate mixture proportions of a mixture g given noisy (error-prone) data from that mixture.
#'
#' @details This is used by the ash function. Most users won't need to call this directly, but is
#' exported for use by some other related packages.
#'
#' @param data list to be passed to log_comp_dens_conv; details depend on model
#' @param g an object representing a mixture distribution (eg normalmix for mixture of normals;
#'  unimix for mixture of uniforms). The component parameters of g (eg the means and variances) specify the
#'  components whose mixture proportions are to be estimated. The mixture proportions of g are the parameters to be estimated;
#'  the values passed in may be used to initialize the optimization (depending on the optmethod used)
#' @param prior numeric vector indicating parameters of "Dirichlet prior"
#'     on mixture proportions
#' @param optmethod name of function to use to do optimization
#' @param control list of control parameters to be passed to optmethod,
#' typically affecting things like convergence tolerance
#' @return list, including the final loglikelihood, the null loglikelihood,
#' an n by k likelihood matrix with (j,k)th element equal to \eqn{f_k(x_j)},
#' the fit
#' and results of optmethod
#' @export
estimate_mixprop = function(data,g,prior,optmethod=c("mixEM","mixVBEM","cxxMixSquarem","mixIP"),control){
  optmethod=match.arg(optmethod)

  pi_init = g$pi
  if(optmethod=="mixVBEM"){pi_init=NULL}  #for some reason pi_init doesn't work with mixVBEM
  k=ncomp(g)

  matrix_llik = t(log_comp_dens_conv(g,data)) #an n by k matrix
  matrix_llik = matrix_llik[!get_exclusions(data),,drop=FALSE] #remove any rows corresponding to excluded cases; saves time in situations where most data are missing
  matrix_llik = matrix_llik - apply(matrix_llik,1, max) #avoid numerical issues by subtracting max of each row
  matrix_lik = exp(matrix_llik)

  # the last of these conditions checks whether the gradient at the null is negative wrt pi0
  # to avoid running the optimization when the global null (pi0=1) is the optimal.
  if(optmethod=="mixVBEM" || max(prior[-1])>1 || min(gradient(matrix_lik)+prior[1]-1,na.rm=TRUE)<0){
    if(optmethod=="cxxMixSquarem"){control=set_control_squarem(control,nrow(matrix_lik))}
    fit=do.call(optmethod,args = list(matrix_lik= matrix_lik, prior=prior, pi_init=pi_init, control=control))
  } else {
    fit = list(converged=TRUE,pihat=c(1,rep(0,k-1)),optmethod="gradient_check")
  }

  ## check if IP method returns negative mixing proportions. If so, run EM.
  if (optmethod == "mixIP" & (min(fit$pihat) < -10 ^ -12)) {
      message("Interior point method returned negative mixing proportions.\n Switching to EM optimization.")
      optmethod <- "mixEM"
      control = list() #use defaults for mixEM in this
      fit = do.call(optmethod, args = list(matrix_lik = matrix_lik,
                                           prior = prior, pi_init = pi_init,
                                           control = control))
  }

  if(!fit$converged){
      warning("Optimization failed to converge. Results may be unreliable. Try increasing maxiter and rerunning.")
  }

  loglik.final =  penloglik(pi,matrix_lik,1) #compute penloglik without penalty
  g$pi=fit$pihat

  return(list(loglik=loglik.final,matrix_lik=matrix_lik,g=g,optreturn=fit,optmethod=optmethod))
}


#' @title Compute Posterior
#'
#' @description Return the posterior on beta given a prior (g) that is
#'     a mixture of normals (class normalmix) and observation
#'     \eqn{betahat ~ N(beta,sebetahat)}
#'
#' @details This can be used to obt
#'
#' @param g a normalmix with components indicating the prior; works
#'     only if g has means 0
#' @param betahat (n vector of observations)
#' @param sebetahat (n vector of standard errors/deviations of
#'     observations)
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
  ldens = stats::dnorm(betahat,0,sqrt(outer(sebetahat^2,sigmaavec^2,FUN="+")),log=TRUE)
  maxldens = apply(ldens, 1, max)
  ldens = ldens - maxldens
  return(exp(ldens))
}

#return the "effective" estimate
#that is the effect size betanew whose z score betanew/se
#would give the same p value as betahat/se compared to a t with df
effective.effect=function(betahat,se,df){
  p = stats::pt(betahat/se,df)
  stats::qnorm(p,sd=se)
}


#' @title Function to compute q values from local false discovery rates
#'
#' @description Computes q values from a vector of local fdr estimates
#'
#' @details The q value for a given lfdr is an estimate of the (tail)
#'     False Discovery Rate for all findings with a smaller lfdr, and
#'     is found by the average of the lfdr for all more significant
#'     findings. See Storey (2003), Annals of Statistics, for
#'     definition of q value.
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
# mode is the location about which inference is going to be centered
# mult is the multiplier by which the sds differ across the grid
autoselect.mixsd = function(data,mult,mode){
  betahat = data$x - mode
  sebetahat = data$s
  exclude = get_exclusions(data)
  betahat = betahat[!exclude]
  sebetahat = sebetahat[!exclude]

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
