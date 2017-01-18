#' @title mean of truncated gamma distribution
#' @description Compute mean of the truncated gamma. 
#' @param a left limit of distribution
#' @param b right limit of distribution
#' @param shape shape of gamma distribution
#' @param rate rate of gamma distribution
#' @export
my_etruncgamma = function(a, b, shape, rate){
  tmp = a
  tmp[a!=b] = (shape/rate*(pgamma(b,shape=shape+1,rate=rate)-pgamma(a,shape=shape+1,rate=rate))/
                 (pgamma(b,shape=shape,rate=rate)-pgamma(a,shape=shape,rate=rate)))[a!=b]
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
  tmp[a!=b] = (shape*(shape+1)/rate^2*(pgamma(b,shape=shape+2,rate=rate)-pgamma(a,shape=shape+2,rate=rate))/
                 (pgamma(b,shape=shape,rate=rate)-pgamma(a,shape=shape,rate=rate)))[a!=b]
  return(tmp)
}