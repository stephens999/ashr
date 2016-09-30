# File contains functions related to logF distribution (cdf,pdf,moments of truncated logF)

#' @title The log-F distribution
#' @description Distribution function for the log-F distribution with \code{df1} and \code{df2}
#'  degrees of freedom (and optional non-centrality parameter \code{ncp}).
#' @param q vector of quantiles
#' @param df1,df2 degrees of freedom
#' @param ncp non-centrality parameter. If omitted the central F is assumed.
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @return The distribution function. 
#' @export
plogf = function(q, df1, df2, ncp, lower.tail=TRUE, log.p=FALSE){
  return(stats::pf(exp(q), df1=df1, df2=df2, ncp=ncp, lower.tail=lower.tail,
            log.p=log.p))
}

#' @title The log-F distribution
#' @description Density function for the log-F distribution with \code{df1} and \code{df2}
#'  degrees of freedom (and optional non-centrality parameter \code{ncp}).
#' @param x vector of quantiles
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @inheritParams plogf
#' @return The density function. 
#' @export
dlogf = function(x, df1, df2, ncp, log=FALSE){
  if (log==FALSE){
    stats::df(exp(x), df1=df1, df2=df2, ncp=ncp)*exp(x)
  }else{
    stats::df(exp(x), df1=df1, df2=df2, ncp=ncp, log=TRUE)+x
  }
}


#' @title my_etrunclogf
#' @description Compute expectation of truncated log-F distribution.
#'
#' @param a Left limit of distribution.
#' @param b Right limit of distribution.
#' @param df1,df2 degrees of freedom
#' @export
my_etrunclogf= function(a,b,df1,df2){
  if (a==b){
    tmp = a
  }else{
    tmp = try(etrunclogf(df1=df1, df2=df2, a=a, b=b, adj=FALSE),silent=TRUE)
    if (class(tmp)=="try-error"){
      tmp = try(etrunclogf(df1=df1, df2=df2, a=a, b=b, adj=TRUE),silent=TRUE)
    }
    
    if (class(tmp)=="try-error"){
      #tmp = NA
      tmp = (a+b)/2
    }
  }
  return(tmp) #deal with extreme case a=b
}

# x*dlogf
etrunclogf_num = function(x,df1,df2,a,b){
  #multiply c to avoid numerical issues
  c = 10^(-round(min(log10(stats::df(exp(a),df1,df2)*exp(a)),
                     log10(stats::df(exp(b),df1,df2)*exp(b)))))
  c*x*stats::df(exp(x),df1=df1,df2=df2)*exp(x)
}

# dlogf
etrunclogf_denom = function(x,df1,df2,a,b){
  #multiply c to avoid numerical issues
  c = 10^(-round(min(log10(stats::df(exp(a),df1,df2)*exp(a)),
                     log10(stats::df(exp(b),df1,df2)*exp(b)))))
  c*stats::df(exp(x),df1=df1,df2=df2)*exp(x)
}

# x multiply by the density of truncated log-F distribution on (a,b) at x
xdtrunclogf = function(x,df1,df2,a,b){
  x*stats::df(exp(x),df1=df1,df2=df2)*exp(x)/(stats::pf(exp(b),df1,df2)-stats::pf(exp(a),df1,df2))
}

# compute expectation of truncated log-F distribution.
etrunclogf = function(df1,df2,a,b,adj=FALSE){
  if (adj==TRUE){
    # numerator and denominator both multiply a constant to avoid numerical issues
    n = stats::integrate(etrunclogf_num, lower=a, upper=b, df1=df1,df2=df2,a=a,b=b)$value
    d = stats::integrate(etrunclogf_denom, lower=a, upper=b, df1=df1,df2=df2,a=a,b=b)$value
    return(n/d)
  }else{
    return(stats::integrate(xdtrunclogf, lower=a, upper=b, df1=df1,df2=df2,a=a,b=b)$value)
  } 
}
