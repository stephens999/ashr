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
  d = dt((x-mean)/sd, df=df, ncp=ncp, log=log)
  if(log){return(d-log(sd))} else{return(d/sd)}
}
