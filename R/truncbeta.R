#' @title mean of truncated Beta distribution
#' @description Compute mean of the truncated Beta. 
#' @param a left limit of distribution
#' @param b right limit of distribution
#' @param alpha,beta shape parameters of Beta distribution
#' @export
my_etruncbeta = function(a, b, alpha, beta){
  tmp = a
  tmp[a!=b] = (alpha/(alpha+beta)*(stats::pbeta(b,alpha+1,beta)-stats::pbeta(a,alpha+1,beta))/
                 (stats::pbeta(b,alpha,beta)-stats::pbeta(a,alpha,beta)))[a!=b]
  # zero denominator case: stats::pbeta(b,alpha,beta) and stats::pbeta(a,alpha,beta) are both 0 or 1 
  tmp[(stats::pbeta(b,alpha,beta)-stats::pbeta(a,alpha,beta))==0] = 
    ifelse(stats::dbeta(a,alpha,beta,log=TRUE)>stats::dbeta(b,alpha,beta,log=TRUE), 
           a, b)[(stats::pbeta(b,alpha,beta)-stats::pbeta(a,alpha,beta))==0]
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
                 (stats::pbeta(b,alpha+2,beta)-stats::pbeta(a,alpha+2,beta))/
                 (stats::pbeta(b,alpha,beta)-stats::pbeta(a,alpha,beta)))[a!=b]
  # zero denominator case: stats::pbeta(b,alpha,beta) and stats::pbeta(a,alpha,beta) are both 0 or 1 
  tmp[(stats::pbeta(b,alpha,beta)-stats::pbeta(a,alpha,beta))==0] = 
    ifelse(stats::dbeta(a,alpha,beta,log=TRUE)>stats::dbeta(b,alpha,beta,log=TRUE), 
           a^2, b^2)[(stats::pbeta(b,alpha,beta)-stats::pbeta(a,alpha,beta))==0]
  return(tmp)
}
