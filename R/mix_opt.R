# This file contains methods for optimizing the mixture proportions given data
# from known mixture components.
# The penalized likelihood being maximized
# is \sum_j \log \sum_k \pi_k f_{jk} + \sum_j (prior_j-1) \log \pi_k

#' @title Estimate mixture proportions of a mixture model by Interior Point method
#'
#' @description Given the individual component likelihoods for a mixture model, estimates the mixture proportions.
#'
#' @details Optimizes \deqn{L(pi)= sum_j w_j log(sum_k pi_k f_{jk}) + h(pi)} 
#' subject to pi_k non-negative and sum_k pi_k = 1. Here \deqn{h(pi)} is
#' a penalty function h(pi) = sum_k (prior_k-1) log pi_k.
#' Calls REBayes::KWDual in the REBayes package, which is in turn a wrapper to the mosek 
#' convex optimization software. So REBayes must be installed to use this. 
#' Used by the ash main function; there is no need for a user to call this 
#' function separately, but it is exported for convenience.
#'
#' 
#' @param matrix_lik, a n by k matrix with (j,k)th element equal to \eqn{f_k(x_j)}.
#' @param prior, a k vector of the parameters of the Dirichlet prior on \eqn{\pi}. Recommended to be rep(1,k)
#' @param pi_init, the initial value of \eqn{\pi} to use. If not specified defaults to (1/k,...,1/k).
#' @param control A list of control parameters to be passed to REBayes::KWDual
#' @param weights weights to be assigned to the observations (an n vector)
#' 
#' @return A list, including the estimates (pihat), the log likelihood for each interation (B)
#' and a flag to indicate convergence
#'  
#' @export
#' 
mixIP = function(matrix_lik, prior, pi_init = NULL, control = list(),
                 weights = NULL) {
 
  # This is the smallest value allowed for the mixture weights.
  min.f <- 0
    
  if(!requireNamespace("REBayes",quietly=TRUE))
    stop("mixIP requires installation of package REBayes")
  control = set_control_mixIP(control)
  n = nrow(matrix_lik)
  k = ncol(matrix_lik)
  
  if(is.null(weights)){weights = rep(1,n)} # give all observations weight 1
  
  A = rbind(diag(length(prior)),matrix_lik) # add in observations corresponding to prior
  w = c(prior-1,weights)
  A = A[w!=0,]    #remove zero weight entries, as these otherwise cause errors
  w = w[w!=0]
  res = REBayes::KWDual(A, rep(1,k), normalize(w), control=control)

  # Fix any mixture weights that are less than the minimum allowed value.
  i <- which(res$f < min.f)
  if (length(i) > 0) {
    warning(paste("Optimization step yields mixture weights that are either",
                  "too small, or negative; weights have been corrected and",
                  "renormalized after the optimization."))
    res$f[i] <- min.f
    res$f    <- normalize(res$f)
  }

  return(list(pihat = normalize(res$f), niter = NULL, converged=(res$status=="OPTIMAL"), control=control))
}

#' @title Estimate mixture proportions of a mixture model using
#' mix-SQP algorithm.
#'
#' @param matrix_lik A matrix containing the conditional likelihood
#' values, possibly normalized.
#'
#' @param prior A vector of the parameters of the Dirichlet prior on
#' the mixture weights.
#'
#' @param pi_init The initial estimate of the mixture weights.
#'
#' @param control A list of settings for the mix-SQP optimization
#' algorithm; see \code{\link[mixsqp]{mixsqp}} for details.
#'
#' @param weights The weights to be assigned to the observations. Must
#' be a vector of length equal the number of rows of \code{matrix_lik}.
#' If \code{weights = NULL}, all observations are assigned the same
#' weight.
#'
#' @return A list object including the estimates (\code{pihat}) and a
#' flag (\code{control}) indicating convergence success or failure.
#' 
#' @importFrom utils modifyList
#' @importFrom mixsqp mixsqp
#' 
#' @export
#' 
mixSQP <- function (matrix_lik, prior, pi_init = NULL,
                    control = list(), weights = NULL) {
  mixsqp.status.converged <- "converged to optimal solution"
  n <- nrow(matrix_lik)
  k <- ncol(matrix_lik)

  # If weights are not provided, set to uniform.
  if (is.null(weights))
    weights <- rep(1,n)

  # If the initial estimate of the mixture weights is not provided,
  # set to uniform.
  if (is.null(pi_init))
    pi_init <- rep(1,k)

  # Add in observations corresponding to the prior.
  A <- rbind(diag(k),matrix_lik) 
  w <- c(prior - 1,weights)
  A <- A[w != 0,,drop=FALSE]
  w <- w[w != 0]

  # Fit the mixture weights using the mix-SQP algorithm.
  control0 <- list(verbose = FALSE)
  control  <- modifyList(control0,control,keep.null = TRUE)
  out      <- mixsqp::mixsqp(A,w,pi_init,control = control)

  # Return the fitted mixture weights, and some other information
  # about the optimization step.
  #
  # Note that the mixture weights may not be normalized (i.e., they
  # may not sum to 1) if the mix-SQP algorithm terminates prematurely,
  # so to be extra cautious, the mix-SQP solution is normalized before
  # returning it.
  return(list(pihat     = out$x/sum(out$x),
              niter     = nrow(out$data),
              converged = (out$status == mixsqp.status.converged),
              control   = control))
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
#' @param matrix_lik, a n by k matrix with (j,k)th element equal to \eqn{f_k(x_j)}.
#' @param prior, a k vector of the parameters of the Dirichlet prior on \eqn{\pi}. Recommended to be rep(1,k)
#' @param pi_init, the initial value of \eqn{\pi} to use. If not specified defaults to (1/k,...,1/k).
#' @param control A list of control parameters for the SQUAREM algorithm, default value is set to be control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE). 
#' 
#' @return A list, including the estimates (pihat), the log likelihood for each interation (B)
#' and a flag to indicate convergence
#'  
#' @export
#' 
#' 
mixEM = function(matrix_lik,prior,pi_init=NULL,control=list()){
  control = set_control_squarem(control,nrow(matrix_lik))
  k=dim(matrix_lik)[2]
  if(is.null(pi_init)){
    pi_init = rep(1/k,k)# Use as starting point for pi
  } 
  res = squarem(par=pi_init,fixptfn=fixpoint, objfn=negpenloglik,matrix_lik=matrix_lik, prior=prior, control=control)
  return(list(pihat = normalize(pmax(0,res$par)), B=res$value.objfn, 
              niter = res$iter, converged=res$convergence, control=control))
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
  return(loglik+penalty(prior, pi))
}

penalty=function(prior, pi){
  subset = (prior != 1.0)
  sum((prior-1)[subset]*log(pi[subset]))
}

#' @title Estimate mixture proportions of a mixture model by EM algorithm (weighted version)
#'
#' @description Given the individual component likelihoods for a mixture model, and a set of weights, estimates the mixture proportions by an EM algorithm.
#'
#' @details Fits a k component mixture model \deqn{f(x|\pi)= \sum_k \pi_k f_k(x)} to independent
#' and identically distributed data \eqn{x_1,\dots,x_n} with weights \eqn{w_1,\dots,w_n}.
#' Estimates mixture proportions \eqn{\pi} by maximum likelihood, or by maximum a posteriori (MAP) estimation for a Dirichlet prior on \eqn{\pi} 
#' (if a prior is specified).  Here the log-likelihood for the weighted data is defined as \eqn{l(\pi) = \sum_j w_j log f(x_j | \pi)}. Uses the SQUAREM package to accelerate convergence of EM. Used by the ash main function; there is no need for a user to call this 
#' function separately, but it is exported for convenience.
#'
#' 
#' @param matrix_lik, a n by k matrix with (j,k)th element equal to \eqn{f_k(x_j)}.
#' @param prior, a k vector of the parameters of the Dirichlet prior on \eqn{\pi}. Recommended to be rep(1,k)
#' @param pi_init, the initial value of \eqn{\pi} to use. If not specified defaults to (1/k,...,1/k).
#' @param weights, an n vector of weights
#' @param control A list of control parameters for the SQUAREM algorithm, default value is set to be control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE). 
#' 
#' @return A list, including the estimates (pihat), the log likelihood for each interation (B)
#' and a flag to indicate convergence
#'  
#' @export
#' 
#' 
w_mixEM = function(matrix_lik,prior,pi_init=NULL,weights=NULL,control=list()){
  control = set_control_squarem(control,nrow(matrix_lik))
  k=dim(matrix_lik)[2]
  if(is.null(pi_init)){
    pi_init = rep(1/k,k)# Use as starting point for pi
  } 
  res = squarem(par=pi_init,fixptfn=w_fixpoint, objfn=w_negpenloglik,matrix_lik=matrix_lik, prior=prior, w=weights,control=control)
  return(list(pihat = normalize(pmax(0,res$par)), B=res$value.objfn, 
              niter = res$iter, converged=res$convergence, control=control))
}

w_fixpoint = function(pi, matrix_lik, prior, w){  
  pi = normalize(pmax(0,pi)) #avoid occasional problems with negative pis due to rounding
  m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  classprob = m/m.rowsum #an n by k matrix
  pinew = normalize(colSums(w*classprob) + prior - 1)
  return(pinew)
}

w_negpenloglik = function(pi,matrix_lik,prior, w){return(-w_penloglik(pi,matrix_lik,prior,w))}

w_penloglik = function(pi, matrix_lik, prior, w){
  pi = normalize(pmax(0,pi))
  m  = t(pi * t(matrix_lik)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  loglik = sum(w*log(m.rowsum))
  subset = (prior != 1.0)
  priordens = sum((prior-1)[subset]*log(pi[subset]))
  return(loglik+priordens)
}


# A vanilla (non-squarem) version of the EM algorithm
# mixEM = function(matrix_lik, prior, pi_init = NULL,tol=0.0001, maxiter=5000){
#   n=nrow(matrix_lik)
#   k=ncol(matrix_lik)
#   B = rep(0,maxiter)
#   pi = pi_init
#   if(is.null(pi_init)){
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
#' @param pi_init the initial value of the posterior parameters. If not specified defaults to the prior parameters.
#' @param control A list of control parameters for the SQUAREM algorithm, default value is set to be   control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE). 
#' 
#' @return A list, whose components include point estimates (pihat), 
#' the parameters of the fitted posterior on \eqn{\pi} (pipost),
#' the bound on the log likelihood for each iteration (B)
#' and a flag to indicate convergence (converged).
#'  
#' @export
#' 
mixVBEM = function(matrix_lik, prior, pi_init = NULL,control=list()){
  control = set_control_squarem(control,nrow(matrix_lik))
  k=ncol(matrix_lik)
  if(is.null(pi_init)){  pi_init = rep(1,k)  }# Use as starting point for pi 
  res = squarem(par=pi_init,fixptfn=VBfixpoint, objfn=VBnegpenloglik,matrix_lik=matrix_lik, prior=prior, control=control)
  
  return(list(pihat = res$par/sum(res$par), B=res$value.objfn, niter = res$iter, converged=res$convergence,post=res$par))
}

# sets up a default for squarem, and modifies it with other provided values
set_control_mixIP=function(control){
  control.default=list(rtol=1e-6)
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  control=utils::modifyList(control.default, control)
  return(control)
}

# sets up a default for squarem, and modifies it with other provided values
set_control_squarem=function(control,nobs){
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  if (nobs > 50000) control.default$trace = TRUE
  control.default$tol = min(0.1/nobs,1.e-7) # set default convergence criteria to be more stringent for larger samples
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  control=utils::modifyList(control.default, control)
  return(control)
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
