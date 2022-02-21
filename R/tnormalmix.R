#' @title Constructor for tnormalmix class
#'
#' @description Creates an object of class tnormalmix (finite mixture
#'   of truncated univariate normals).
#'
#' @param pi Cector of mixture proportions (length k say).
#' 
#' @param mean Vector of means (length k).
#' 
#' @param sd Vector of standard deviations (length k).
#' 
#' @param a Vector of left truncation points of each component (length k).
#' 
#' @param b Cector of right truncation points of each component (length k).
#'
#' @return An object of class \dQuote{tnormalmix}.
#'
#' @export
#'
#' @examples tnormalmix(c(0.5,0.5),c(0,0),c(1,2),c(-10,0),c(0,10))
#'
tnormalmix = function (pi, mean, sd, a, b)
  structure(data.frame(pi,mean,sd,a,b),class = "tnormalmix")

#' @title comp_sd.normalmix
#' 
#' @description Returns standard deviations of the truncated normal mixture.
#' 
#' @param m A truncated normal mixture distribution with k components.
#' 
#' @return A vector of length k.
#' 
#' @export
#' 
comp_sd.tnormalmix = function(m)
  sqrt(my_vtruncnorm(m$a,m$b,m$mean,m$sd))

#' @title comp_mean.tnormalmix
#' 
#' @description Returns mean of the truncated-normal mixture.
#' 
#' @param m A truncated normal mixture distribution with k components.
#' 
#' @return A vector of length k.
#' 
#' @export
#' 
comp_mean.tnormalmix = function (m)
  my_etruncnorm(m$a,m$b,m$mean,m$sd)

#' @importFrom stats dnorm
#' @importFrom stats pnorm
#' 
#' @export
#' 
comp_dens.tnormalmix = function (m, y, log = FALSE) {
  k = ncomp(m)
  n = length(y)
  d = matrix(rep(y,rep(k,n)),nrow = k)
  # No cases of b = a.
  return(matrix(dnorm(d,m$mean,m$sd))/(pnorm(m$b) - pnorm(m$a)))
}

#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' 
#' @export
#' 
comp_dens_conv.tnormalmix = function (m, data) {
  if (!is_normal(data$lik))
    stop("Error: truncated normal mixture for non-normal likelihood is not ",
         "yet implemented")
  if (length(data$s) == 1)
    data$s = rep(data$s,length(data$x))
  A      = sqrt(outer(1/m$sd^2,1/data$s^2,FUN = "+"))
  B      = 1/sqrt(outer(m$sd^2,data$s^2,FUN = "+"))
  C      = outer(m$sd,data$s,"/")
  D      = pnorm(m$b/m$sd) - pnorm(m$a/m$sd)
  varmat = outer(m$sd^2,data$s^2,FUN = "+")
  left   = pnorm(A*m$b - t(t(B*C)*data$x))
  right  = pnorm(A*m$a - t(t(B*C)*data$x))
  denx   = dnorm(matrix(data$x,length(m$sd),length(data$x),
                        byrow = TRUE)/sqrt(varmat))/sqrt(varmat)
  result = ((left - right)*denx)/D
  DD     = dnorm(m$b/m$sd)
  lleft  = dnorm(A*m$b - t(t(B*C)*data$x))
  result[m$a == m$b,] = ((lleft*denx/varmat)/DD)[m$a == m$b,]
  result[m$sd == 0,] = denx[m$sd == 0,]
  return(result)
}

#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' 
#' @export
#' 
log_comp_dens_conv.tnormalmix = function (m, data) {
  if (!is_normal(data$lik))
    stop("Error: truncated normal mixture for non-normal likelihood is not ",
         "yet implemented")
  
  # Use previous function directly.
  if (length(data$s) == 1)
    data$s = rep(data$s,length(data$x))
  A = sqrt(outer(1/m$sd^2,1/data$s^2,FUN = "+"))
  B = 1/sqrt(outer(m$sd^2,data$s^2,FUN = "+"))
  C = outer(m$sd,data$s,"/")
  D = pnorm(m$b/m$sd) - pnorm(m$a/m$sd)
  varmat = outer(m$sd^2,data$s^2,FUN = "+")
  left   = pnorm(A*m$b - t(t(B*C)*data$x))
  right  = pnorm(A*m$a - t(t(B*C)*data$x))
  denx   = dnorm(t(matrix(data$x,length(data$x),length(m$sd))),0,sqrt(varmat),
                 log = TRUE)
  result = log(left - right) + denx -
           log(matrix(D,length(m$sd),length(data$x)))
  DD     = dnorm(m$b/m$sd)
  lleft  = dnorm(A*m$b - t(t(B*C)*data$x))
  result[m$a  == m$b,] = (log(lleft/DD) + denx - log(varmat))[m$a == m$b,]
  result[m$sd == 0,]   = denx[m$sd == 0,]
  return(result)
}

#' @importFrom stats pnorm
#' 
#' @export
#' 
comp_cdf.tnormalmix = function (m, y, lower.tail = TRUE) {
  k = length(m$pi)
  n = length(y)
  tmp = matrix(1,nrow = k,ncol = n)
  subset = outer(m$a,y,">")
  tmp[subset] = 0
  subset1 = outer(m$a,y,"<=")
  subset2 = outer(m$b,y,">=")
  subset  = subset1 & subset2
  if (sum(subset) > 0) {
    YY  = matrix(y,k,n,byrow = TRUE)
    MM  = matrix(m$mean,k,n)
    SD  = matrix(m$sd,k,n)
    pnc = matrix(pnorm(YY[subset],MM[subset],SD[subset]),k,n)
    A   = matrix(m$a,k,ncol = n)
    pna = matrix(pnorm(A[subset],MM[subset],SD[subset],lower.tail),k,n)
    B   = matrix(m$b,k,ncol = n)
    pnb = matrix(pnorm(B[subset],MM[subset],SD[subset],lower.tail),k,n)
  }
  tmp[subset] = (pnc-pna)/(pnb-pna)
  return(tmp)
}

#' @importFrom stats pnorm
#' 
#' @export
#' 
comp_cdf_post.tnormalmix = function (m, c, data) {
  if (!is_normal(data$lik)) 
    stop("Error: truncated normal mixture for non-normal likelihood is ",
         "not yet implemented")
  k = length(m$pi)
  n = length(data$x)
  tmp = matrix(1,nrow = k,ncol = n)
  tmp[m$a > c,] = 0
  subset = m$a <= c & m$b >= c
  if (sum(subset)>0) {
    X   = 1/(outer(data$s^2,m$sd[subset]^2,FUN = "/") + 1)
    Y   = 1/outer(1/data$s^2,1/m$sd[subset]^2,FUN = "+")
    A   = matrix(m$a[subset],nrow = sum(subset),ncol = n)
    pna = pnorm(t(A),X*data$x + t(t(1-X) * m$mean[subset]),sqrt(Y))
    C   = matrix(c,nrow = sum(subset),ncol = n)
    pnc = pnorm(t(C),X*data$x + t(t(1-X) * m$mean[subset]),sqrt(Y))
    B   = matrix(m$b[subset],nrow = sum(subset),ncol = n)
    pnb = pnorm(t(B),X*data$x + t(t(1-X) * m$mean[subset]),sqrt(Y))
  }
  tmp[subset,] = t((pnc - pna)/(pnb - pna))
  subset       = (m$a == m$b)
  tmp[subset,] = rep(m$a[subset] <= c,n)
  subset       = B == C
  tmp[subset]  = 1
  ### ZMZ: in case of pnc = pnb, we make it 1 and other
  ### Nan 0 to eliminate the 0/0.
  ### use naive situation
  tmp[is.nan(tmp)] = 0
  return(tmp)
}

#' @export
comp_postmean.tnormalmix = function (m, data) {
  if (!is_normal(data$lik))
    stop("Error: truncated normal mixture for non-normal likelihood is not ",
         "yet implemented")
  k = length(m$pi)
  n = length(data$x)
  A = 1/(outer(m$sd^2,data$s^2,FUN = "/") + 1)
  B = 1/outer(1/m$sd^2,1/data$s^2,FUN = "+")
  result = my_etruncnorm(m$a,m$b,A*m$mean + t(t(1-A)*data$x),sqrt(B))
  ismissing = which(is.na(data$x) | is.na(data$s))
  if (length(ismissing) > 0)
    result[,ismissing] = m$mean
  return(result)
}

#' @export
comp_postsd.tnormalmix = function (m, data) {
  if (!is_normal(data$lik))
    stop("Error: truncated normal mixture for non-normal likelihood is not ",
         "yet implemented")
  k = length(m$pi)
  n = length(data$x)
  A = 1/(outer(m$sd^2,data$s^2,FUN = "/") + 1)
  B = 1/outer(1/m$sd^2,1/data$s^2,FUN = "+")
  result = sqrt(my_vtruncnorm(m$a,m$b,t(A*m$mean + t(t(1 - A)*data$x)),
                              t(sqrt(B))))
  return(t(result))
}

#' @export
comp_postmean2.tnormalmix = function (m, data)
  comp_postsd.tnormalmix(m,data)^2 + comp_postmean.tnormalmix(m,data)^2

