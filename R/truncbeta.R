#' @title mean of truncated Beta distribution
#' @description Compute mean of the truncated Beta. 
#' @param a left limit of distribution
#' @param b right limit of distribution
#' @param alpha,beta shape parameters of Beta distribution
#' @export
my_etruncbeta = function(a, b, alpha, beta){
  tmp = a
  tmp[a!=b] = (alpha/(alpha+beta)*(pbeta(b,alpha+1,beta)-pbeta(a,alpha+1,beta))/
                 (pbeta(b,alpha,beta)-pbeta(a,alpha,beta)))[a!=b]
  # zero denominator case: pbeta(b,alpha,beta) and pbeta(a,alpha,beta) are both 0 or 1 
  tmp[(pbeta(b,alpha,beta)-pbeta(a,alpha,beta))==0] = 
    ifelse(dbeta(a,alpha,beta,log=TRUE)>dbeta(b,alpha,beta,log=TRUE), 
           a, b)[(pbeta(b,alpha,beta)-pbeta(a,alpha,beta))==0]
  return(tmp)
}

#' @title second moment of truncated Beta distribution
#' @description Compute second moment of the truncated Beta. 
#' @param a left limit of distribution
#' @param b right limit of distribution
#' @param alpha,beta shape parameters of Beta distribution
#' @export
my_e2truncbeta = function(a, b, alpha, beta){
  tmp = a^2
  tmp[a!=b] = (alpha*(alpha+1)/((alpha+beta)*(alpha+beta+1))*
                 (pbeta(b,alpha+2,beta)-pbeta(a,alpha+2,beta))/
                 (pbeta(b,alpha,beta)-pbeta(a,alpha,beta)))[a!=b]
  # zero denominator case: pbeta(b,alpha,beta) and pbeta(a,alpha,beta) are both 0 or 1 
  tmp[(pbeta(b,alpha,beta)-pbeta(a,alpha,beta))==0] = 
    ifelse(dbeta(a,alpha,beta,log=TRUE)>dbeta(b,alpha,beta,log=TRUE), 
           a^2, b^2)[(pbeta(b,alpha,beta)-pbeta(a,alpha,beta))==0]
  return(tmp)
}