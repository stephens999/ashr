#' @title mean of truncated gamma distribution
#' @description Compute mean of the truncated gamma. 
#' @param a left limit of distribution
#' @param b right limit of distribution
#' @param shape shape of gamma distribution
#' @param rate rate of gamma distribution
#' @export
my_etruncgamma = function(a, b, shape, rate){
  tmp = a
  tmp[a!=b] = (shape/rate*(stats::pgamma(b,shape=shape+1,rate=rate)-stats::pgamma(a,shape=shape+1,rate=rate))/
                 (stats::pgamma(b,shape=shape,rate=rate)-stats::pgamma(a,shape=shape,rate=rate)))[a!=b]
  # zero denominator case: stats::pgamma(b,shape,rate) and stats::pgamma(a,shape,rate) are both 0 or 1 
  tmp[(stats::pgamma(b,shape=shape,rate=rate)-stats::pgamma(a,shape=shape,rate=rate))==0] = 
    ifelse(stats::dgamma(a,shape=shape,rate=rate,log=TRUE)>stats::dgamma(b,shape=shape,rate=rate,log=TRUE), 
           a, b)[(stats::pgamma(b,shape=shape,rate=rate)-stats::pgamma(a,shape=shape,rate=rate))==0]
  return(tmp)
}

#' @title second moment of truncated gamma distribution
#' @description Compute second moment of the truncated gamma. 
#' @param a left limit of distribution
#' @param b right limit of distribution
#' @param shape shape of gamma distribution
#' @param rate rate of gamma distribution
#' @export
my_e2truncgamma = function(a, b, shape, rate){
  tmp = a^2
  tmp[a!=b] = (shape*(shape+1)/rate^2*(stats::pgamma(b,shape=shape+2,rate=rate)-stats::pgamma(a,shape=shape+2,rate=rate))/
                 (stats::pgamma(b,shape=shape,rate=rate)-stats::pgamma(a,shape=shape,rate=rate)))[a!=b]
  # zero denominator case: stats::pgamma(b,shape,rate) and stats::pgamma(a,shape,rate) are both 0 or 1 
  tmp[(stats::pgamma(b,shape=shape,rate=rate)-stats::pgamma(a,shape=shape,rate=rate))==0] = 
    ifelse(stats::dgamma(a,shape=shape,rate=rate,log=TRUE)>stats::dgamma(b,shape=shape,rate=rate,log=TRUE), 
           a^2, b^2)[(stats::pgamma(b,shape=shape,rate=rate)-stats::pgamma(a,shape=shape,rate=rate))==0]
  return(tmp)
}
