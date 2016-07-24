#' @title Return lfsr, lfdr, etc from ash object
#'
#' @description These functions simply return elements of an ash object, generally without doing any calculations. 
#' (So if the value was not computed during the original call to ash, eg because of how outputlevel was set in the call, 
#' then NULL will be returned.)
#' Accessing elements in this way
#' rather than directly from the ash object will help ensure compatability moving forward
#' (e.g. if the internal structure of the ash object changes during software development.)
#'
#' @param a the fitted ash object
#'
#' @describeIn get_lfsr local false sign rate
#' @export
get_lfsr=function(a){a$result$lfsr}

#' @describeIn get_lfsr local false discovery rate
#' @export
get_lfdr=function(a){a$result$lfdr}

#' @describeIn get_lfsr svalue
#' @export
get_svalue=function(a){a$result$svalue}

#' @describeIn get_lfsr qvalue
#' @export
get_qvalue=function(a){a$result$qvalue}

#' @describeIn get_lfsr posterior mean
#' @export
get_pm=function(a){a$result$PosteriorMean}

#' @describeIn get_lfsr posterior standard deviation
#' @export
get_psd=function(a){a$result$PosteriorSD}

#' @describeIn get_lfsr positive probability
#' @export
get_pp=function(a){a$result$PositiveProb}

#' @describeIn get_lfsr negative probability
#' @export
get_np=function(a){a$result$NegativeProb}

#' @describeIn get_lfsr log-likelihood 
#' @export
get_loglik=function(a){a$loglik}

#' @describeIn get_lfsr log-likelihood ratio
#' @export
get_logLR=function(a){a$logLR}

#' @describeIn get_lfsr fitted g mixture
#' @export
get_fitted_g=function(a){a$fitted_g}


#' @describeIn get_lfsr pi0, the proportion of nulls 
#' @export
get_pi0 = function(a){
  null.comp = comp_sd(a$fitted_g)==0
  return(sum(a$fitted_g$pi[null.comp]))
}
