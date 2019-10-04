#' @useDynLib ashr
#' @import Matrix truncnorm SQUAREM doParallel pscl Rcpp foreach parallel
#' @title Adaptive Shrinkage
#'
#' @description Implements Empirical Bayes shrinkage and false discovery rate
#' methods based on unimodal prior distributions.
#'
#' @details The ash function provides a number of ways to perform Empirical Bayes shrinkage
#' estimation and false discovery rate estimation. The main assumption is that
#' the underlying distribution of effects is unimodal. Novice users are recommended
#' to start with the examples provided below.
#'
#' In the simplest case the inputs to ash are a vector of estimates (betahat)
#' and their corresponding standard errors (sebetahat), and degrees of freedom (df).
#' The method assumes that for some (unknown) "true" vector of effects beta, the statistic
#' (betahat[j]-beta[j])/sebetahat[j] has a $t$ distribution on $df$ degrees of freedom.
#' (The default of df=NULL assumes a normal distribution instead of a t.)
#'
#' By default the method estimates the vector beta under the assumption that beta ~ g for a distribution
#' g in G, where G is some unimodal family of distributions to be specified (see parameter \code{mixcompdist}).
#' By default is to assume the mode is 0, and this is suitable for settings where you are interested in testing which beta[j]
#' are non-zero. To estimate the mode see parameter \code{mode}.
#'
#' As is standard in empirical Bayes methods, the fitting proceeds in two stages:
#' i) estimate g by maximizing a (possibly penalized) likelihood;
#' ii) compute the posterior distribution for each beta[j] | betahat[j],sebetahat[j]
#' using the estimated g as the prior distribution.
#'
#' A more general case allows that beta[j]/sebetahat[j]^alpha | sebetahat[j] ~ g.
#'
#' @param betahat a p vector of estimates
#'
#' @param sebetahat a p vector of corresponding standard errors
#'
#' @param mixcompdist distribution of components in mixture used to represent the family G.
#' Depending on the choice of mixture component, the family G becomes more or less flexible.
#' Options are:\cr
#' \describe{
#' \item{uniform}{G is (approximately) any symmetric unimodal distribution}
#' \item{normal}{G is (approximately) any scale mixture of normals}
#' \item{halfuniform}{G is (approximately) any unimodal distribution}
#' \item{+uniform}{G is (approximately) any unimodal distribution with support constrained to be greater than the mode.}
#' \item{-uniform}{G is (approximately) any unimodal distribution with support constrained to be less than the mode.}
#' \item{halfnormal}{G is (approximately) any scale mixture of truncated normals where the normals are truncated at the mode}
#' }
#' If you are happy to assume a symmetric distribution for effects, you can use
#' "uniform" or "normal". If you believe your effects
#' may be asymmetric, use "halfuniform" or "halfnormal". If you want
#' to allow only positive/negative effects use "+uniform"/"-uniform".
#' The use of "normal" and "halfnormal" is permitted only if df=NULL.
#'
#' @param df appropriate degrees of freedom for (t) distribution of
#' (betahat-beta)/sebetahat; default is NULL which is actually treated as
#' infinity (Gaussian)
#'
#' @param method specifies how ash is to be run. Can be "shrinkage"
#' (if main aim is shrinkage) or "fdr" (if main aim is to assess false discovery rate
#' or false sign rate (fsr)). This is simply a convenient way to specify certain
#' combinations of parameters: "shrinkage" sets pointmass=FALSE and
#' prior="uniform"; "fdr" sets pointmass=TRUE and prior="nullbiased".
#'
#' @param optmethod specifies the function implementing an
#' optimization method.
#'
#' @param nullweight scalar, the weight put on the prior under
#' "nullbiased" specification, see \code{prior}
#'
#' @param mode either numeric (indicating mode of g) or string
#' "estimate", to indicate mode should be estimated, or a two
#' dimension numeric vector to indicate the interval to be searched
#' for the mode.
#'
#' @param pointmass Logical, indicating whether to use a point mass at
#' zero as one of components for a mixture distribution.
#'
#' @param prior string, or numeric vector indicating Dirichlet prior
#' on mixture proportions (defaults to "uniform", or (1,1...,1); also
#' can be "nullbiased" (nullweight,1,...,1) to put more weight on
#' first component), or "unit" (1/K,...,1/K) [for optmethod=mixVBEM
#' version only].
#'
#' @param mixsd Vector of standard deviations for underlying mixture components.
#'
#' @param gridmult the multiplier by which the default grid values for
#' mixsd differ by one another. (Smaller values produce finer grids.)
#'
#' @param outputlevel Determines amount of output. There are several
#' numeric options: 0 = just fitted g; 1 = also PosteriorMean and
#' PosteriorSD; 2 = everything usually needed; 3 = also include results
#' of mixture fitting procedure (including matrix of log-likelihoods
#' used to fit mixture). 4 and 5 are reserved for outputting additional
#' data required by the (in-development) flashr package. The user can
#' also specify the output they require in detail (see Examples).
#'
#' @param g The prior distribution for beta. Usually this is unspecified (NULL) and
#' estimated from the data. However, it can be used in conjuction with fixg=TRUE
#' to specify the g to use (e.g. useful in simulations to do computations with the "true" g).
#' Or, if g is specified but fixg=FALSE, the g specifies the initial value of g used before optimization,
#' (which also implicitly specifies mixcompdist).
#'
#' @param fixg If TRUE, don't estimate g but use the specified g -
#' useful for computations under the "true" g in simulations.
#'
#' @param alpha Numeric value of alpha parameter in the model.
#'
#' @param grange Two dimension numeric vector indicating the left and
#' right limit of g. Default is c(-Inf, Inf).
#'
#' @param control A list of control parameters passed to optmethod.
#'
#' @param lik Contains details of the likelihood used; for general
#' ash. Currently, the following choices are allowed: normal (see
#' function lik_normal(); binomial likelihood (see function
#' lik_binom); likelihood based on logF error distribution (see
#' function lik_logF); mixture of normals likelihood (see function
#' lik_normalmix); and Poisson likelihood (see function lik_pois).
#'
#' @param weights a vector of weights for observations; use with
#' optmethod = "w_mixEM"; this is currently beta-functionality.
#'
#' @param pi_thresh a threshold below which to prune out mixture
#' components before computing summaries (speeds up computation since
#' empirically many components are usually assigned negligible
#' weight). The current implementation still returns the full fitted
#' distribution; this only affects the posterior summaries.
#'
#' @param ... Further arguments of function \code{ash} to be passed to
#' \code{\link{ash.workhorse}}.
#'
#' @return ash returns an object of \code{\link[base]{class}} "ash", a
#' list with some or all of the following elements (determined by
#' outputlevel) \cr
#' \item{fitted_g}{fitted mixture}
#' \item{loglik}{log P(D|fitted_g)}
#' \item{logLR}{log[P(D|fitted_g)/P(D|beta==0)]}
#' \item{result}{A dataframe whose columns are:}
#' \describe{
#'   \item{NegativeProb}{A vector of posterior probability that beta is
#'     negative.}
#'   \item{PositiveProb}{A vector of posterior probability that beta is
#'     positive.}
#'   \item{lfsr}{A vector of estimated local false sign rate.}
#'   \item{lfdr}{A vector of estimated local false discovery rate.}
#'   \item{qvalue}{A vector of q values.}
#'   \item{svalue}{A vector of s values.}
#'   \item{PosteriorMean}{A vector consisting the posterior mean of beta
#'     from the mixture.}
#'   \item{PosteriorSD}{A vector consisting the corresponding posterior
#'     standard deviation.}
#'   }
#' \item{call}{a call in which all of the specified arguments are
#'   specified by their full names}
#' \item{data}{a list containing details of the data and models
#'   used (mostly for internal use)}
#' \item{fit_details}{a list containing results of mixture optimization,
#'   and matrix of component log-likelihoods used in this optimization}
#'
#' @seealso \code{\link{ashci}} for computation of credible intervals
#' after getting the ash object return by \code{ash()}
#'
#' @export ash
#' @export ash.workhorse
#'
#' @examples
#'
#' beta = c(rep(0,100),rnorm(100))
#' sebetahat = abs(rnorm(200,0,1))
#' betahat = rnorm(200,beta,sebetahat)
#' beta.ash = ash(betahat, sebetahat)
#' names(beta.ash)
#' head(beta.ash$result) # the main dataframe of results
#' head(get_pm(beta.ash)) # get_pm returns posterior mean
#' head(get_lfsr(beta.ash)) # get_lfsr returns the local false sign rate
#' graphics::plot(betahat,get_pm(beta.ash),xlim=c(-4,4),ylim=c(-4,4))
#'
#' \dontrun{
#' # Why is this example included here? -Peter
#' CIMatrix=ashci(beta.ash,level=0.95)
#' print(CIMatrix)
#' }
#'
#' # Illustrating the non-zero mode feature.
#' betahat=betahat+5
#' beta.ash = ash(betahat, sebetahat)
#' graphics::plot(betahat,get_pm(beta.ash))
#' betan.ash=ash(betahat, sebetahat,mode=5)
#' graphics::plot(betahat,get_pm(betan.ash))
#' summary(betan.ash)
#'
#' # Running ash with different error models
#' beta.ash1 = ash(betahat, sebetahat, lik = lik_normal())
#' beta.ash2 = ash(betahat, sebetahat, lik = lik_t(df=4))
#'
#' e = rnorm(100)+log(rf(100,df1=10,df2=10)) # simulated data with log(F) error
#' e.ash = ash(e,1,lik=lik_logF(df1=10,df2=10))
#'
#' # Specifying the output
#' beta.ash = ash(betahat, sebetahat, output = c("fitted_g","logLR","lfsr"))
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
#'
#' # running with weights
#' beta.ash = ash(betahat, sebetahat, optmethod="w_mixEM",
#'                weights = c(rep(0.5,100),rep(1,100)))
#'
#' # Different algorithms can be used to compute maximum-likelihood
#' # estimates of the mixture weights. Here, we illustrate use of the
#' # EM algorithm and the (default) SQP algorithm.
#' set.seed(1)
#' betahat  <- c(8.115,9.027,9.289,10.097,9.463)
#' sebeta   <- c(0.6157,0.4129,0.3197,0.3920,0.5496)
#' fit.em   <- ash(betahat,sebeta,mixcompdist = "normal",optmethod = "mixEM")
#' fit.sqp  <- ash(betahat,sebeta,mixcompdist = "normal",optmethod = "mixSQP")
#' range(fit.em$fitted$pi - fit.sqp$fitted$pi)
ash <- function (betahat, sebetahat,
                 mixcompdist = c("uniform","halfuniform","normal","+uniform",
                                 "-uniform","halfnormal"),
                 df = NULL,...){
  # This calls ash.workhorse, but then modifies the returned list so that the call is the original ash call
  a = ash.workhorse(betahat,sebetahat,mixcompdist = mixcompdist,df = df,...)
  utils::modifyList(a,list(call = match.call()))
}

#' @describeIn ash Adaptive Shrinkage with full set of options.
ash.workhorse <-
    function(betahat, sebetahat, method = c("fdr","shrink"),
             mixcompdist = c("uniform","halfuniform","normal","+uniform",
                             "-uniform","halfnormal"),
             optmethod = c("mixSQP","mixIP","cxxMixSquarem","mixEM",
                           "mixVBEM","w_mixEM"),
             df = NULL,nullweight = 10,pointmass = TRUE,
             prior = c("nullbiased","uniform","unit"),mixsd = NULL,
             gridmult = sqrt(2),outputlevel = 2,g = NULL,fixg = FALSE,
             mode = 0,alpha = 0,grange = c(-Inf,Inf),control = list(),lik = NULL, weights=NULL, pi_thresh = 1e-10) {

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
  if(!is.numeric(prior)){
    prior = match.arg(prior)
  }

  ## Check to see whether df is Inf. If so, switch to NULL.
  if (length(df) > 1) {
    stop("Only one value can be specified for df.")
  }
  if (!is.null(df) && is.infinite(df)) {
    df <- NULL
  }

  # set likelihood based on defaults if missing
  if(is.null(lik)){
    if(is.null(df)){
      lik = lik_normal()
    } else {lik = lik_t(df)}
  }

  # poisson likelihood has non-negative g
  # do not put big weight on null component
  # automatically estimate the mode if not specified
  if(lik$name=="pois"){
    if (lik$data$link=="identity"){
      grange = c(max(0,min(grange)), max(grange))
    }
    if(missing(nullweight)){nullweight = 1}
    if(missing(mode) & missing(g)){mode = "estimate"}
  }
  # binomial likelihood (identity link) has g restricted on [0,1]
  if(lik$name=="binom"){
    if (lik$data$link=="identity"){
      grange = c(max(0,min(grange)), min(1,max(grange)))
    }
    if(missing(nullweight)){nullweight = 1}
    if(missing(mode) & missing(g)){mode = "estimate"}
  }

  #if length of mode is 2 then the numeric values give the range of values to search
  if(sum(mode=="estimate") | length(mode)==2){ #just pass everything through to ash.estmode for non-zero-mode
    args <- as.list(environment())
    args$mode = NULL
    args$outputlevel = NULL
    args$method=NULL # avoid specifying method as well as prior/pointmass
    args$g = NULL # avoid specifying g as well as mode
    mode = ifelse(is.numeric(mode),mode,NA)

    # set range to search the mode
    if (lik$name=="pois"){
      if (lik$data$link=="identity"){
        args$modemin = min(mode, min(lik$data$y/lik$data$scale),na.rm = TRUE)
        args$modemax = max(mode, max(lik$data$y/lik$data$scale),na.rm = TRUE)
      }else if (lik$data$link=="log"){
        args$modemin = min(log(lik$data$y/lik$data$scale+0.01))
        args$modemax = max(log(lik$data$y/lik$data$scale+0.01))
      }
    }else if(lik$name=="binom"){
      if (lik$data$link=="identity"){
        args$modemin = min(grange)
        args$modemax = max(grange)
      }else if (lik$data$link=="logit"){
        logitp = log((lik$data$y+0.01)/(lik$data$n+0.02)/(1-(lik$data$y+0.01)/(lik$data$n+0.02)))
        args$modemin = min(logitp)
        args$modemax = max(logitp)
      }

    }else{
      args$modemin = min(mode, min(betahat),na.rm = TRUE)
      args$modemax = max(mode, max(betahat),na.rm = TRUE)
    }
    mode = do.call(ash.estmode,args)
  }

  ##1.Handling Input Parameters
  mixcompdist = match.arg(mixcompdist)
  optmethod   = match.arg(optmethod)

  # Set optimization method
  optmethod = set_optmethod(optmethod)
  check_args(mixcompdist,df,prior,optmethod,gridmult,sebetahat,betahat)
  check_lik(lik, betahat, sebetahat, df, mixcompdist) # minimal check that it obeys requirements
  lik = add_etruncFUN(lik) #if missing, add a function to compute mean of truncated distribution
  data = set_data(betahat, sebetahat, lik, alpha)

  ##2. Generating mixture distribution g

  if (fixg & missing(g)) {
    stop("If fixg = TRUE then you must specify g!")
  }

  
  if (!is.null(g)) {
    k = ncomp(g)
    null.comp = 1 # null.comp not actually used
    prior = setprior(prior, k, nullweight, null.comp)
  } else {
    if (mixcompdist %in% c("uniform", "halfuniform", "+uniform", "-uniform")) {
      # For unimix prior, if mode is exactly the boundary of g's range, have
      #   to use "+uniform" or "-uniform"
      if (min(grange) == mode) {
        mixcompdist = "+uniform"
      } else if (max(grange) == mode) {
        mixcompdist = "-uniform"
      }
    }

    if (is.null(mixsd)) {
      mixsd = autoselect.mixsd(data, gridmult, mode, grange, mixcompdist)
    }
    if (pointmass) {
      mixsd = c(0, mixsd)
    }
    null.comp = which.min(mixsd) # which component is the "null"

    k = length(mixsd)
    prior = setprior(prior, k, nullweight, null.comp)
    pi = initpi(k, length(data$x), null.comp)

    if (mixcompdist == "normal")
      g = normalmix(pi, rep(mode, k), mixsd)
    if (mixcompdist == "uniform")
      g = unimix(pi, mode - mixsd, mode + mixsd)
    if (mixcompdist == "+uniform")
      g = unimix(pi, rep(mode, k), mode + mixsd)
    if (mixcompdist == "-uniform")
      g = unimix(pi, mode - mixsd, rep(mode, k))
    if (mixcompdist == "halfuniform") {
      if (min(mixsd) > 0) { #no point mass
        # Simply reflect the components.
        g = unimix(c(pi, pi) / 2,
                   c(mode - mixsd, rep(mode, k)),
                   c(rep(mode, k), mode + mixsd))
        prior = c(prior, prior)
      } else {
        # Define two sets of components, but don't duplicate null component.
        null.comp = which.min(mixsd)
        g = unimix(c(pi, pi[-null.comp]) / (2 - pi[null.comp]),
                   c(mode - mixsd, rep(mode, k-1)),
                   c(rep(mode, k), (mode + mixsd)[-null.comp]))
        prior = c(prior, prior[-null.comp])
      }
    }
    if (mixcompdist == "halfnormal") {
      if (min(mixsd) > 0) { #no point mass
        g = tnormalmix(c(pi, pi) / 2,
                       rep(mode, 2 * k),
                       c(mixsd, mixsd),
                       c(rep(-Inf, k), rep(0, k)),
                       c(rep(0, k), rep(Inf, k)))
        prior = c(prior, prior)
      } else {
        null.comp = which.min(mixsd)
        g = tnormalmix(c(pi, pi[-null.comp]) / (2 - pi[null.comp]),
                       rep(mode, 2 * k - 1),
                       c(mixsd, mixsd[-null.comp]),
                       c(rep(-Inf, k), rep(0, k - 1)),
                       c(rep(0, k), rep(Inf, k - 1)))
        prior = c(prior, prior[-null.comp])
      }
    }

    # Constrain g within grange.
    gconstrain = constrain_mix(g, prior, grange, mixcompdist)
    g = gconstrain$g
    prior = gconstrain$prior
  }

  #check that all prior are >=1 (as otherwise have problems with infinite penalty)
  if(!all(prior>=1) & optmethod != "mixVBEM"){
    stop("Error: prior must all be >=1 (unless using optmethod mixVBEM)")}

  
  ##3. Fitting the mixture
  if(!fixg){
    pi.fit=estimate_mixprop(data,g,prior,optmethod=optmethod,control=control,weights=weights)
  } else {
    pi.fit = list(g=g,penloglik = calc_loglik(g,data)+penalty(prior, g$pi))
  }

  ##4. Computing the return values

  val = list() # val will hold the return value
  ghat = pi.fit$g
  output = set_output(outputlevel) #sets up flags for what to output
  if("flash_data" %in% output){
    flash_data = calc_flash_data(ghat, data, pi.fit$penloglik)
    val = c(val, list(flash_data = flash_data))
  }
  if("fitted_g" %in% output){val = c(val,list(fitted_g=ghat))}
  if("loglik" %in% output){val = c(val,list(loglik =calc_loglik(ghat,data)))}
  if("logLR" %in% output){val = c(val,list(logLR=calc_logLR(ghat,data)))}
  if("data" %in% output){val = c(val,list(data=data))}
  if("fit_details" %in% output){val = c(val,list(fit_details = pi.fit))}
  if("post_sampler" %in% output){
    val = c(val,list(post_sampler=function(nsamp){post_sample(ghat, data, nsamp)}))
  }

  # Compute the result component of value -
  # result is a dataframe containing lfsr, etc
  # resfns is a list of functions used to produce columns of that dataframe
  resfns = set_resfns(output)
  if(length(resfns)>0){
    result = data.frame(betahat = betahat,sebetahat = sebetahat)
    if(!is.null(df)){result$df = df}
    result = cbind(result,as.data.frame(lapply(resfns,do.call,list(g=prune(ghat,pi_thresh),data=data))))
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
#' @export
compute_lfsr = function(NegativeProb,ZeroProb){
  ifelse(NegativeProb> 0.5*(1-ZeroProb),1-NegativeProb,NegativeProb+ZeroProb)
}



#The kth element of this vector is the derivative
#of the loglik for $\pi=(\pi_0,...,1-\pi_0,...)$ with respect to $\pi_0$ at $\pi_0=1$.
gradient = function(matrix_lik){
  n = nrow(matrix_lik)
  grad = n - ColsumModified(matrix_lik)
  return(grad)
}

ColsumModified = function(matrix_l){
  small = abs(matrix_l) < 10e-100
  matrix_l[small] = matrix_l[small]+10e-100
  colSums(matrix_l/matrix_l[,1])
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
#' @param weights vector of weights (for use with w_mixEM; in beta)
#' @return list, including the final loglikelihood, the null loglikelihood,
#' an n by k likelihood matrix with (j,k)th element equal to \eqn{f_k(x_j)},
#' the fit
#' and results of optmethod
#'
#' @export
#'
estimate_mixprop = function (data, g, prior,
                             optmethod = c("mixSQP", "mixEM", "mixVBEM",
                                           "cxxMixSquarem", "mixIP", "w_mixEM"),
                             control, weights = NULL) {

  optmethod = match.arg(optmethod)

  pi_init = g$pi
  # In rare cases, a mixture proportion that is initialized to zero can cause
  #   optimization to fail. So we ensure that all proportions are positive.
  if (!all(pi_init > 0)) {
    pi_init = pmax(pi_init, 1e-6)
    pi_init = pi_init / sum(pi_init)
  }
  # For some reason pi_init doesn't work with mixVBEM.
  if (optmethod=="mixVBEM") {
    pi_init = NULL
  }

  matrix_llik = t(log_comp_dens_conv(g, data)) # an n by k matrix
  # Remove excluded cases; saves time when most data is missing.
  matrix_llik = matrix_llik[!get_exclusions(data), , drop=FALSE]
  # Avoid numerical issues by subtracting the max of each row.
  lnorm = apply(matrix_llik, 1, max)
  matrix_llik = matrix_llik - lnorm
  matrix_lik = exp(matrix_llik)

  # All-zero columns pose problems for most optimization methods.
  nonzero_cols = (apply(matrix_lik, 2, max) > 0)
  if (!all(nonzero_cols)) {
    prior = prior[nonzero_cols]
    weights = weights[nonzero_cols]
    pi_init = pi_init[nonzero_cols]
    matrix_lik = matrix_lik[, nonzero_cols, drop=FALSE]
  }

  ncomponents = length(prior)

  if (!is.null(weights) && !(optmethod %in% c("w_mixEM", "mixIP", "mixSQP")))
    stop("Weights can only be used with optmethod w_mixEM, mixIP or mixSQP.")
  if (optmethod %in% c("w_mixEM", "mixSQP") && is.null(weights)) {
    weights = rep(1, nrow(matrix_lik))
  }
  if (ncomponents > 1 && !is.null(weights)) {
    fit = do.call(optmethod, args = list(matrix_lik = matrix_lik,
                                         prior = prior,
                                         pi_init = pi_init,
                                         control = control,
                                         weights = weights))
  } else if (ncomponents > 1
             && (optmethod == "mixVBEM"
                 || max(prior[-1]) > 1
                 || min(gradient(matrix_lik) + prior[1] - 1, na.rm = TRUE) < 0)) {
    # The last condition checks whether the gradient at the null is negative
    #   wrt pi0. This avoids running the optimization when the global null
    #   (pi0 = 1) is optimal.
    if (optmethod == "cxxMixSquarem") {
      control = set_control_squarem(control, nrow(matrix_lik))
    }
    fit = do.call(optmethod, args = list(matrix_lik = matrix_lik,
                                         prior = prior,
                                         pi_init = pi_init,
                                         control = control))
  } else {
    fit = list(converged = TRUE,
               pihat = c(1, rep(0, ncomponents - 1)),
               optmethod = "gradient_check")
  }

  if (!fit$converged) {
      warning("Optimization failed to converge. Results may be unreliable. ",
              "Try increasing maxiter and rerunning.")
  }

  fit$pihat = pmax(fit$pihat, 0)

  g$pi = rep(0, ncomp(g))
  g$pi[nonzero_cols] = fit$pihat
  # Value of objective function:
  penloglik = penloglik(fit$pihat, matrix_lik, prior) + sum(lnorm)

  return(list(penloglik = penloglik,
              matrix_lik = matrix_lik,
              g = g,
              optreturn = fit,
              optmethod = optmethod))
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
# grange is the user-specified range of mixsd
autoselect.mixsd = function(data,mult,mode,grange,mixcompdist){
  if (data$lik$name %in% c("pois")){
    data$x = data$lik$data$y/data$lik$data$scale #estimate of lambda
    data$s = sqrt(data$x)/data$lik$data$scale #standard error of estimate
    # if the link is log we probably want to take the log of this?
  }
  if (data$lik$name %in% c("binom")){
    data$x = data$lik$data$y
  }
  
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

# Constrain g within grange.
# g: unimix, normalmix, or tnormalmix prior
# prior: k-vector
# grange: two-dimensional numeric vector indicating the left and right limit
#   of the prior g.
constrain_mix = function(g, prior, grange, mixcompdist) {
  pi = g$pi
  if (inherits(g, "normalmix") || inherits(g, "tnormalmix")) {
    # Normal mixture priors always lie on (-Inf, Inf), so ignore grange.
    if (max(grange) < Inf | min(grange) > -Inf) {
      warning("Can't constrain grange for normal/halfnormal mixture prior ",
              "case.")
    }
  } else if (inherits(g, "unimix")) {
    # Truncate the uniform mixture components that are out of grange.
    g$a = pmax(g$a, min(grange))
    g$b = pmin(g$b, max(grange))
    compidx = !duplicated(cbind(g$a, g$b)) # remove duplicated components
    pi = pi[compidx]
    if (sum(pi) == 0) {
      stop("No component has positive mixture probability after constraining ",
           "range.")
    }
    pi = pi / sum(pi)
    g = unimix(pi, g$a[compidx], g$b[compidx])
    prior = prior[compidx]
  } else {
    stop("constrain_mix does not recognize that prior type.")
  }

  return(list(g = g, prior = prior))
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
