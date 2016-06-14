# This file contains methods for optimizing the mixture proportions given data
# from known mixture components.
# The penalized likelihood being maximized
# is \sum_j \log \sum_k \pi_k f_{jk} + \sum_j (prior_j-1) \log \pi_k

#' @title Estimate mixture proportions of a mixture model by Interior Point method
#'
#' @description Given the individual component likelihoods for a mixture model, estimates the mixture proportions.
#'
#' @details Fits a k component mixture model \deqn{f(x|\pi)= \sum_k \pi_k f_k(x)} to independent
#' and identically distributed data \eqn{x_1,\dots,x_n}. 
#' Estimates mixture proportions \eqn{\pi} by maximum likelihood, or by maximum a posteriori (MAP) estimation for a Dirichlet prior on \eqn{\pi} 
#' (if a prior is specified). Calls REBayes::KWDual in the REBayes package, which is in turn a wrapper to the mosek 
#' convex optimization software. So REBayes must be installed to use this. Used by the ash main function; there is no need for a user to call this 
#' function separately, but it is exported for convenience.
#'
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
mixIP = function(matrix_lik, prior, pi_init = NULL, control = list()){
  if(!requireNamespace("REBayes",quietly=TRUE)){stop("mixIP requires installation of package REBayes")}
  n = nrow(matrix_lik)
  k = ncol(matrix_lik)
  #A = matrix_lik
  A = rbind(diag(length(prior)),matrix_lik) # add in observations corresponding to prior
  w = c(prior-1,rep(1,n))
  A = A[w!=0,]    #remove zero weight entries, as these otherwise cause errors
  w = w[w!=0]
  #w = rep(1,n+k)
  res = REBayes::KWDual(A, rep(1,k), normalize(w), control=control)
  return(list(pihat = normalize(res$f), niter = NULL, converged=(res$status=="OPTIMAL")))
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
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)
  
  k=dim(matrix_lik)[2]
  if(is.null(pi_init)){
    pi_init = rep(1/k,k)# Use as starting point for pi
  } 
  res = squarem(par=pi_init,fixptfn=fixpoint, objfn=negpenloglik,matrix_lik=matrix_lik, prior=prior, control=controlinput)
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
  control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
  namc=names(control)
  if (!all(namc %in% names(control.default))) 
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput=modifyList(control.default, control)
  
  k=ncol(matrix_lik)
  if(is.null(pi_init)){  pi_init = rep(1,k)  }# Use as starting point for pi 
  res = squarem(par=pi_init,fixptfn=VBfixpoint, objfn=VBnegpenloglik,matrix_lik=matrix_lik, prior=prior, control=controlinput)
  
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

