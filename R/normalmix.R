############################### METHODS FOR normalmix class ###########################

#' @title Constructor for normalmix class
#'
#' @description Creates an object of class normalmix (finite mixture
#'     of univariate normals)
#'
#' @details None
#'
#' @param pi vector of mixture proportions
#' @param mean vector of means
#' @param sd vector of standard deviations
#'
#' @return an object of class normalmix
#'
#' @export
#'
#' @examples normalmix(c(0.5,0.5),c(0,0),c(1,2))
#'
normalmix = function(pi,mean,sd){
  structure(data.frame(pi,mean,sd),class="normalmix")
}


#' @title comp_sd.normalmix
#' @description returns sds of the normal mixture
#' @param m a normal mixture distribution with k components
#' @return a vector of length k
#' @export
comp_sd.normalmix = function(m){
  m$sd
}

#' @title comp_mean.normalmix
#' @description returns mean of the normal mixture
#' @param m a normal mixture distribution with k components
#' @return a vector of length k
#' @export
comp_mean.normalmix = function(m){
  m$mean
}

#' @export
comp_dens.normalmix = function(m,y,log=FALSE){
  k=ncomp(m)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(stats::dnorm(d, m$mean, m$sd, log),nrow=k))
}

#' @export
comp_cdf.normalmix = function(m,y,lower.tail=TRUE){
  vapply(y,stats::pnorm,m$mean,m$mean,m$sd,lower.tail)
}



#' @title comp_dens_conv.normalmix
#' @description returns density of convolution of each component of a
#'     normal mixture with N(0,s^2) at x. Note that
#'     convolution of two normals is normal, so it works that way
#' @param m mixture distribution with k components
#' @param data a list with components x and s to be interpreted as a normally-distributed observation and its standard error
#' @return a k by n matrix
comp_dens_conv.normalmix = function(m,data){
  if(!is_normal(data$lik)){
    stop("Error: normal mixture for non-normal likelihood is not yet implemented")
  }
  sdmat = sqrt(outer(data$s^2,m$sd^2,FUN="+")) #n by k matrix of standard deviations of convolutions
  return(t(stats::dnorm(outer(data$x,m$mean,FUN="-")/sdmat)/sdmat))
}

#' @title log_comp_dens_conv.normalmix
#' @description returns log-density of convolution of each component
#'     of a normal mixture with N(0,s^2) or s*t(v) at x. Note that
#'     convolution of two normals is normal, so it works that way
#' @inheritParams comp_dens_conv.normalmix
#' @return a k by n matrix
log_comp_dens_conv.normalmix = function(m,data){
  if(!is_normal(data$lik)){
    stop("Error: normal mixture for non-normal likelihood is not yet implemented")
  }
  sdmat = sqrt(outer(data$s^2,m$sd^2,"+")) #n by k matrix of standard deviations of convolutions
  return(t(stats::dnorm(outer(data$x,m$mean,FUN="-")/sdmat,log=TRUE) - log(sdmat)))
}


#' @title comp_cdf_conv.normalmix
#' @description returns cdf of convolution of each component of a
#'     normal mixture with N(0,s^2) at x. Note that
#'     convolution of two normals is normal, so it works that way
#' @param m mixture distribution with k components
#' @param data a list with components x and s to be interpreted as a normally-distributed observation and its standard error
#' @return a k by n matrix
comp_cdf_conv.normalmix = function (m, data) {
  if(!is_normal(data$lik)){
    stop("Error: normal mixture for non-normal likelihood is not yet implemented")
  }
  sdmat = sqrt(outer(data$s^2, m$sd^2, FUN="+")) #n by k matrix of standard deviations of convolutions
  return(t(stats::pnorm(outer(data$x, m$mean, FUN="-") / sdmat)))
}


#' @export
comp_cdf_post.normalmix=function(m,c,data){
  if(!is_normal(data$lik)){
    stop("Error: normal mixture for non-normal likelihood is not yet implemented")
  }
  k = length(m$pi)
  
  #compute posterior standard deviation (s1) and posterior mean (m1)
  s1 = sqrt(outer(data$s^2,m$sd^2,FUN="*")/outer(data$s^2,m$sd^2,FUN="+"))
  ismissing = (is.na(data$x) | is.na(data$s))
  s1[ismissing,]=m$sd
  
  m1 = t(comp_postmean(m,data))
  t(stats::pnorm(c,mean=m1,sd=s1))
}


#' @export
comp_postmean.normalmix = function(m,data){
  if(!is_normal(data$lik)){
    stop("Error: normal mixture for non-normal likelihood is not yet implemented")
  }
  tmp=(outer(data$s^2,m$mean, FUN="*") + outer(data$x,m$sd^2, FUN="*"))/outer(data$s^2,m$sd^2,FUN="+")
  ismissing = (is.na(data$x) | is.na(data$s))
  tmp[ismissing,]=m$mean #return prior mean when missing data
  t(tmp)
}

#' @export
comp_postsd.normalmix = function(m,data){
  if(!is_normal(data$lik)){
    stop("Error: normal mixture for non-normal likelihood is not yet implemented")
  }
  t(sqrt(outer(data$s^2,m$sd^2,FUN="*")/outer(data$s^2,m$sd^2,FUN="+")))
}

#' @export
comp_postmean2.normalmix = function(m,data){
  comp_postsd(m,data)^2 + comp_postmean(m,data)^2
}

#' @title post_sample.normalmix
#' 
#' @description returns random samples from the posterior, given a
#'   prior distribution m and n observed datapoints.
#' 
#' @param m mixture distribution with k components
#' @param data a list with components x and s to be interpreted as a 
#'     normally-distributed observation and its standard error
#' @param nsamp number of samples to return for each observation
#' @return a nsamp by n matrix
#' @importFrom stats rnorm
#' @export
post_sample.normalmix = function(m,data,nsamp){
  k = length(m$pi)
  n = length(data$x)
  
  postprob = comp_postprob(m,data)
  postmean = comp_postmean(m,data)
  postsd = comp_postsd(m,data)
  
  # Sample mixture components
  mixcomp = apply(postprob, 2, function(prob) {
    sample(1:k, nsamp, replace=TRUE, prob=prob)
  })
  # Use samples to index into postmean and postsd matrices
  idx = mixcomp + rep(k*(0:(n-1)), each=nsamp)
  samp = rnorm(nsamp*n, postmean[idx], postsd[idx])
  matrix(samp, nrow=nsamp, ncol=n)
}
