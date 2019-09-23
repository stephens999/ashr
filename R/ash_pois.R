#' @title Performs adaptive shrinkage on Poisson data
#' @description Uses Empirical Bayes to fit the model \deqn{y_j | \lambda_j ~ Poi(c_j \lambda_j)} with \deqn{h(lambda_j) ~ g()}
#' where \eqn{h} is a specified link function (either "identity" or "log" are permitted). 

#' @details The model is fit in two stages: i) estimate \eqn{g} by maximum likelihood (over the set of symmetric
#' unimodal distributions) to give estimate \eqn{\hat{g}}; 
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{y_j,\hat{g}}.
#'    Note that the link function \eqn{h} affects the prior assumptions (because, e.g., assuming a unimodal prior on \eqn{\lambda} is
#'    different from assuming unimodal on \eqn{\log\lambda}), but posterior quantities are always computed for the
#'    for \eqn{\lambda} and *not* \eqn{h(\lambda)}.
#' @param y vector of Poisson observations.
#' @param scale vector of scale factors for Poisson observations: the model is \eqn{y[j]~Pois(scale[j]*lambda[j])}.
#' @param link string, either "identity" or "log", indicating the link function. 
#' @param ... other parameters to be passed to ash
#' 
#' @examples 
#'    beta = c(rep(0,50),rexp(50))
#'    y = rpois(100,beta) # simulate Poisson observations
#'    y.ash = ash_pois(y,scale=1)
#' @export
ash_pois = function(y, scale=1, link=c("identity","log"), ...){
  link = match.arg(link)
  ash(rep(0,length(y)), 1, lik=lik_pois(y,scale,link), ...)
}