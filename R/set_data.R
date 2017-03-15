# Issues: 
# - the cdfFUN has to take log parameter currently (used in log_comp_dens_conv)
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
#' @param lik a likelihood (see e.g., lik_normal())
#' @param alpha specifies value of alpha to use (model is for betahat/sebetahat^alpha | sebetahat)
#' 
#' @return data object (list) 
#' @export
set_data = function(betahat, sebetahat, lik=NULL, alpha=0){
 
  if(length(sebetahat)==1L){sebetahat = rep(sebetahat, length(betahat))}
  
  data=list(x = betahat/(sebetahat^alpha),
            s = sebetahat^(1-alpha),
            alpha=alpha,
            s_orig = sebetahat)

  if(is.null(lik)){lik = lik_normal()}
  data$lik = lik

  return(data)
}

#extract data corresponding to ith data point
extract_data=function(data,i){
  if(!is_const(data$lik)){stop("extracting data not supported for non-constant likelihoods")}
  data_i = list(x=data$x[i],
                s=data$s[i],
                s_orig = data$s_orig[i],
                alpha = data$alpha,
                lik = data$lik)
  return(data_i)
}

n_obs = function(data){return(length(data$x))}

get_exclusions=function(data){
  return((data$s==0 | data$s == Inf | is.na(data$x) | is.na(data$s)))
}
