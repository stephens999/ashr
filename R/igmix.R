
############################### METHODS FOR igmix class ###########################

#' @title Constructor for igmix class
#'
#' @description Creates an object of class igmix (finite mixture of
#'     univariate inverse-gammas)
#'
#' @details None
#'
#' @param pi vector of mixture proportions
#' @param alpha vector of shape parameters
#' @param beta vector of rate parameters
#'
#' @return an object of class igmix
#'
#' @export
#'
#' @examples igmix(c(0.5,0.5),c(1,1),c(1,2))
#'
igmix = function(pi,alpha,beta){
  structure(data.frame(pi,alpha,beta),class="igmix")
}

comp_sd.igmix = function(m){
  m$beta/(m$alpha-1)/sqrt(m$alpha-2)
}

#' @export
comp_mean.igmix = function(m){
  m$beta/(m$alpha-1)
}

comp_dens.igmix = function(m,y,log=FALSE){
  k=ncomp(m)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(stats::dgamma(1/d, shape=m$alpha, rate=outer(m$beta,1/y^2),log),nrow=k))
}

#density of product of each component of a inverse-gamma mixture with Gamma(v/2,v/2) at s
# s an n-vector at which density is to be evaluated
#return a k by n matrix
comp_dens_conv.igmix = function(m,data,FUN="+"){
  k=ncomp(m)
  x = data$x
  s = data$s
  v = data$v
  n=length(s)
  dens = t(exp(v/2*log(v/2)-lgamma(v/2)
               +(v/2-1)*outer(log(s^2),rep(1,k))
               +outer(rep(1,n),m$alpha*log(m$beta)-lgamma(m$alpha)+lgamma(m$alpha+v/2))
               -outer(rep(1,n),m$alpha+v/2)*log(outer(v/2*s^2,m$beta,FUN="+"))))
  return(dens)
  
}

#' @export
comp_cdf.igmix = function(m,y,lower.tail=TRUE){
  vapply(y,pigamma,m$alpha,m$alpha,m$beta,lower.tail)
}


#' @export
comp_cdf_post.igmix=function(m,c,data){
  #compute posterior shape (alpha1) and rate (beta1)
  alpha1 = m$alpha+data$v/2
  beta1 = outer(m$beta,data$v/2*data$s^2,FUN="+")
  ismissing = is.na(data$s)
  beta1[,ismissing]=m$beta
  return(t(pigamma(c,alpha1,beta1)))
}


#' @export
comp_postmean.igmix = function(m,data){
  k = length(m$pi)
  n=length(data$s)
  tmp=outer(data$v/2*data$s^2,m$beta,FUN="+")/outer(rep(1,n),m$alpha+data$v/2-1)
  ismissing = is.na(data$s)
  tmp[ismissing,]=m$beta/(m$alpha-1) #return prior mean when missing data
  t(tmp)
}


#' @export
comp_postsd.igmix = function(m,data){
  k = length(m$pi)
  n=length(data$s)
  alpha1 = outer(rep(1,n),m$alpha+data$v/2-1)
  beta1 = outer(data$v/2*data$s^2,m$beta,FUN="+")
  return(t(beta1/(alpha1-1)/sqrt(alpha1-2)))
}

#' @export
comp_dens.igmix = function(m,y,log=FALSE){
  k=ncomp(m)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  ig_dens=matrix(densigamma(d, m$alpha, m$beta),nrow=k)
  if(log==TRUE){ig_dens=log(ig_dens)}
  return(ig_dens)
}
