

############################### METHODS FOR unimix class ###########################

#' @title Constructor for unimix class
#'
#' @description Creates an object of class unimix (finite mixture of
#'     univariate uniforms)
#'
#' @details None
#'
#' @param pi vector of mixture proportions
#' @param a vector of left hand ends of uniforms
#' @param b vector of right hand ends of uniforms
#'
#' @return an object of class unimix
#'
#' @export
#'
#' @examples unimix(c(0.5,0.5),c(0,0),c(1,2))
unimix = function(pi,a,b){
  structure(data.frame(pi,a,b),class="unimix")
}

#' @export
comp_cdf.unimix = function(m,y,lower.tail=TRUE){
  vapply(y,stats::punif,m$a,min=m$a,max=m$b,lower.tail)
}

comp_sd.unimix = function(m){
  (m$b-m$a)/sqrt(12)
}

#' @export
comp_mean.unimix = function(m){
  (m$a+m$b)/2
}



comp_dens.unimix = function(m,y,log=FALSE){
  k=ncomp(m)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(stats::dunif(d, m$a, m$b, log),nrow=k))
}

#' density of convolution of each component of a unif mixture 
#' @param m a mixture of class unimix
#' @param data, see set_data()
#'
#' @return a k by n matrix
#'
#' @export
comp_dens_conv.unimix = function(m,data){
  return(exp(log_comp_dens_conv(m,data)))
}

#' log density of convolution of each component of a unif mixture 
#' @inheritParams comp_dens_conv.unimix
#' @return a k by n matrix of densities
log_comp_dens_conv.unimix = function(m,data){
  b = pmax(m$b,m$a) #ensure a<b
  a = pmin(m$b,m$a)
  lik = data$lik
  
  lpa = do.call(lik$lcdfFUN, list(outer(data$x,a,FUN="-")/data$s))
  lpb = do.call(lik$lcdfFUN, list(outer(data$x,b,FUN="-")/data$s))
  
  if (sum(lpa-lpb,na.rm=TRUE)<0){
    tmp = lpa
    lpa = lpb
    lpb = tmp
  }
   
  lcomp_dens = t(lpa + log(1-exp(lpb-lpa))) - log(b-a)
  lcomp_dens[a==b,] = t(do.call(lik$lpdfFUN, list(outer(data$x,b,FUN="-")/data$s))
                       -log(data$s))[a==b,]
  return(lcomp_dens)
}

#' cdf of convolution of each component of a unif mixture 
#' @param m a mixture of class unimix
#' @param data, see set_data()
#'
#' @return a k by n matrix
#'
#' @export
#'
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' 
comp_cdf_conv.unimix = function (m, data) {
  if(!is_normal(data$lik)){
    stop("Error: diagnostic plot of uniform mixture for non-normal likelihood is not yet implemented")
  }
  b = pmax(m$b,m$a) #ensure a<b
  a = pmin(m$b,m$a)
  b.mat = outer(data$x, b, FUN = "-") / data$s
  a.mat = outer(data$x, a, FUN = "-") / data$s
  lcomp_cdf_conv = t(log(pmax(a.mat * pnorm(a.mat) + dnorm(a.mat) - b.mat * pnorm(b.mat) - dnorm(b.mat), 0)) - log(a.mat - b.mat))
  lcomp_cdf_conv[a == b, ] = t(log(pnorm(a.mat)))[a == b, ]
  return(exp(lcomp_cdf_conv))
}


#' @export
comp_cdf_post.unimix=function(m,c,data){
  k = length(m$pi)
  n=length(data$x)
  lik = data$lik
  
  tmp = matrix(1,nrow=k,ncol=n)
  tmp[m$a > c,] = 0
  subset = m$a<=c & m$b>c # subset of components (1..k) with nontrivial cdf
  if(sum(subset)>0){
    lpna = do.call(lik$lcdfFUN, list(outer(data$x,m$a[subset],FUN="-")/data$s))
    lpnc = do.call(lik$lcdfFUN, list(outer(data$x,rep(c,sum(subset)),FUN="-")/data$s))
    lpnb = do.call(lik$lcdfFUN, list(outer(data$x,m$b[subset],FUN="-")/data$s))
    tmp[subset,] = t((exp(lpnc-lpna)-1)/(exp(lpnb-lpna)-1))
    #tmp[subset,] = t((pnc-pna)/(pnb-pna)) ; doing on different log scale reduces numerical issues
  }
  subset = (m$a == m$b) #subset of components with trivial cdf
  tmp[subset,]= rep(m$a[subset] <= c,n)
  #Occasionally we would encounter issue such that in some entries pna[i,j]=pnb[i,j]=pnc[i,j]=0 or pna=pnb=pnc=1
  #Those are the observations with significant betahat(small sebetahat), resulting in pnorm() return 1 or 0
  #due to the thin tail property of normal distribution.(or t-distribution, although less likely to occur)
  #Then R would be dividing 0 by 0, resulting  in NA values
  #In practice, those observations would have 0 probability of belonging to those "problematic" components
  #Thus any sensible value in [0,1] would not matter much, as they are highly unlikely to come from those
  #components in posterior distribution.
  #Here we simply assign the "naive" value as as (c-a)/(b-a)
  #As the component pdf is rather smaller over the region.
  tmpnaive=matrix(rep((c-m$a)/(m$b-m$a),length(data$x)),nrow=k,ncol=n)
  tmp[is.nan(tmp)]= tmpnaive[is.nan(tmp)]
  tmp
}

#note that with uniform prior, posterior is truncated normal, so
#this is computed using formula for mean of truncated normal
#' @export
comp_postmean.unimix = function(m,data){
  x=data$x
  s=data$s

  lik = data$lik
  
  alpha = outer(x, -m$b,FUN="+")/s
  beta = outer(x, -m$a, FUN="+")/s
  
  tmp = x-s*do.call(lik$etruncFUN, list(alpha,beta))
  
  # alpha = outer(-x, m$a,FUN="+")/s
  # beta = outer(-x, m$b, FUN="+")/s
  # 
  # tmp = x + s*do.call(lik$etruncFUN, list(alpha,beta))
 
  ismissing = is.na(x) | is.na(s)
  tmp[ismissing,]= (m$a+m$b)/2
  t(tmp)
}

# as for posterior mean, but compute posterior mean squared value
#' @export
comp_postmean2.unimix = function(m,data){
  x=data$x
  s=data$s

  lik = data$lik
  alpha = outer(-x, m$a,FUN="+")/s
  beta = outer(-x, m$b, FUN="+")/s
  tmp = x^2 +
    2*x*s* do.call(lik$etruncFUN, list(alpha,beta)) +
    s^2* do.call(lik$e2truncFUN, list(alpha,beta)) 
 
  ismissing = is.na(x) | is.na(s)
  tmp[ismissing,]= (m$b^2+m$a*m$b+m$a^2)/3
  t(tmp)
}

# #not yet implemented!
# #just returns 0s for now
# comp_postsd.unimix = function(m,data){
#   k= ncomp(m)
#   n=length(data$x)
#   return(matrix(NA,nrow=k,ncol=n))
#   #  return(sqrt(comp_postmean2(m,betahat,sebetahat,v)-comp_postmean(m,betahat,sebetahat,v)^2))
# }


#' @title post_sample.unimix
#' 
#' @description returns random samples from the posterior, given a
#'   prior distribution m and n observed datapoints.
#' 
#' @param m mixture distribution with k components
#' @param data a list with components x and s to be interpreted as a 
#'     normally-distributed observation and its standard error
#' @param nsamp number of samples to return for each observation
#' @return a nsamp by n matrix
#' @importFrom truncnorm rtruncnorm
#' @export
post_sample.unimix = function(m,data,nsamp){
  k = length(m$pi)
  n = length(data$x)
  
  postprob = comp_postprob(m,data)
  # Sample mixture components
  mixcomp = apply(postprob, 2, function(prob) {
    sample(1:k, nsamp, replace=TRUE, prob=prob)
  })
  
  a = m$a[mixcomp]
  b = m$b[mixcomp]
  
  samp = rtruncnorm(nsamp*n, a = a, b = b, 
                    mean = rep(data$x, each=nsamp), 
                    sd = rep(data$s, each=nsamp))
  # rtruncnorm gives NA when a = b, so these need to be set separately:
  idx = (a == b)
  samp[idx] = a[idx]
  
  matrix(samp, nrow=nsamp, ncol=n)
}
