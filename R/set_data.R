# Issues: 
# - the cdfFUN has to take log parameter currently (used in log_compdens_conv)
# - the FUNargs have to be the same to cdf, pdf and etruncFUN; this means I added mean and sd to my_etrunct, but that functino
#   actually doesn't use them. Need to think about what this means.
# - 

#' Takes raw data and sets up data object for use by ash
#' 
#' @details The data object stores both the data, and details of the model to be used for the data.
#' For example, in the generalized version of ash the cdf and pdf of the likelihood are
#' stored here.
#' 
#' @param betahat vector of betahats
#' @param sebetahat vector of standard errors
#' @param df degree of freedom to assume (for t likelihood)
#' @param lik a likelihood (see eg normal_lik())
#' @param alpha specifies value of alpha to use (model is for betahat/sebetahat^alpha | sebetahat)
#' 
#' @return data object (list) 
#' @export
set_data = function(betahat, sebetahat, lik=NULL, alpha=0){
  # Dealing with precise input of betahat, currently we exclude them from the EM algorithm
  if(length(sebetahat)==1L){sebetahat = rep(sebetahat, length(betahat))}
  exclude = get_exclusions(betahat,sebetahat)
  #exclude = rep(FALSE,length(betahat))
  
  data=list()
  data$x = betahat[!exclude]/(sebetahat[!exclude]^alpha)
  data$s = sebetahat[!exclude]^(1-alpha)
  
  data$exclude = exclude
  data$alpha=alpha
  data$s_orig = sebetahat[!exclude]
  
  if(is.null(lik)){lik = normal_lik()}
  data$lik = lik

  return(data)
}

get_exclusions=function(betahat,sebetahat){
  return((sebetahat==0 | sebetahat == Inf | is.na(betahat) | is.na(sebetahat)))
}