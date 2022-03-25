#' @title Expected Value of Truncated Normal
#' 
#' @description Computes the means of truncated normal distributions with
#'   parameters \code{a}, \code{b}, \code{mean}, and \code{sd}. Arguments
#'   can be scalars, vectors, or matrices. Arguments of shorter length will
#'   be recycled according to the usual recycling rules, but \code{a} and 
#'   \code{b} must have the same length. Missing values are accepted for all 
#'   arguments.
#'
#' @param a The lower limit for the support of the truncated normal. Can be
#'   \code{-Inf}.
#' 
#' @param b The upper limit for the support. Can be \code{Inf}. \code{a} and 
#'   \code{b} must have the same length, and each element of \code{a} should 
#'   be less than or equal to the corresponding element of \code{b}.
#' 
#' @param mean The mean of the untruncated normal.
#' 
#' @param sd The standard deviation of the untruncated normal. Standard
#'   deviations of zero are interpreted as numerically (rather than exactly)
#'   zero, so that the untruncated mean is returned if it lies within 
#'   \code{[a, b]} and the nearer of \code{a} and \code{b} is returned
#'   otherwise.
#' 
#' @return The expected values of truncated normal distributions with
#'   parameters \code{a}, \code{b}, \code{mean}, and \code{sd}. If any of the
#'   arguments is a matrix, then a matrix will be returned.
#'   
#' @seealso \code{\link{my_e2truncnorm}}, \code{\link{my_vtruncnorm}}
#'   
#' @export
my_etruncnorm = function(a, b, mean = 0, sd = 1) {
  # based off of https://github.com/cossio/TruncatedNormal.jl
  do_truncnorm_argchecks(a, b)
  
  # initialize array to store result
  if (is.matrix(a)){
    res = array(dim = dim(0 * (a + b + mean + sd)))
  }
  else{
    res = rep(NA, length.out = length(0 * (a + b + mean + sd)))
  }
    
  #make sure input sizes all match 
  a = rep(a, length.out = length(res))
  b = rep(b, length.out = length(res))
  mean = rep(mean, length.out = length(res))
  sd = rep(sd, length.out = length(res))
  
  # Handle zero sds. Return the mean of the untruncated normal when it is 
  #   located inside of the interval [alpha, beta]. Otherwise, return the 
  #   endpoint that is closer to the mean.
  sd.zero = (sd == 0)
  res[sd.zero & b <= mean] = b[sd.zero & b <= mean]
  res[sd.zero & a >= mean] = a[sd.zero & a >= mean]
  res[sd.zero & a < mean & b > mean] = mean[sd.zero & a < mean & b > mean]
  
  # Handle NAN inputs
  isna = is.na(a) | is.na(b) | is.na(mean) | is.na(sd)
  res[isna] = NA
  
  # Focus in on where sd is nonzero and nothing is nan
  a = a[!sd.zero & !isna]
  b = b[!sd.zero & !isna]
  mean = mean[!sd.zero & !isna]
  sd = sd[!sd.zero & !isna]
  
  # Rescale to standard normal distributions
  alpha = (a - mean) / sd
  beta = (b - mean) / sd
  # initialize array for scaled 2nd moments
  scaled.mean = rep(NA, length.out = length(alpha))
  
  # point mass bc endpoints are equal
  # 1st moment is point-value
  endpoints_equal = (alpha == beta)
  scaled.mean[endpoints_equal] = alpha[endpoints_equal]
  # keep track of which spots in scaled.2mom are already computed
  computed = endpoints_equal
  
  # force to satisfy β ≥ 0 and |α| ≤ |β|
  # so either α ≤ 0 ≤ β or 0 < α ≤ β
  flip = !computed & abs(alpha) > abs(beta)
  flip[is.na(flip)] = FALSE
  orig.alpha = alpha
  alpha[flip] = -beta[flip]
  beta[flip] = -orig.alpha[flip]
  
  # both endpoints infinite/untruncated normal distribution
  # scaled mean is 0
  both_inf = !computed & is.infinite(alpha) & is.infinite(beta)
  scaled.mean[both_inf] = 0
  computed = computed | both_inf
  
  # truncated to [α,∞) 
  # 2nd moment simplifies to ϕ(α)/(1 - Φ(α))
  beta_inf = !computed & is.infinite(beta)
  scaled.mean[beta_inf] = sqrt(2/pi) / Re(erfcx(alpha[beta_inf] / sqrt(2)))
  computed = computed | beta_inf
  
  # a ≤ 0 ≤ b
  #catestrophic cancellation is less of an issue
  alpha_negative = !computed & alpha <= 0
  diff = (beta[alpha_negative] - alpha[alpha_negative]) * (alpha[alpha_negative] + beta[alpha_negative]) / 2
  #√(2/π) * expm1(-Δ) * exp(-α^2 / 2) / erf(β/√2, α/√2)
  scaled.mean[alpha_negative] = sqrt(2/pi) * expm1(-diff) * exp(-alpha[alpha_negative]^2 / 2) 
  denom = Re(erf(alpha[alpha_negative] / sqrt(2))) - Re(erf(beta[alpha_negative] / sqrt(2)))
  scaled.mean[alpha_negative] = scaled.mean[alpha_negative] / denom
  computed = computed | alpha_negative
  
  # 0 < a < b
  #strategically avoid catestrophic cancellation as much as possible
  diff = (beta[!computed] - alpha[!computed]) * (alpha[!computed] + beta[!computed]) / 2
  denom = exp(-diff) * Re(erfcx(beta[!computed] / sqrt(2))) - Re(erfcx(alpha[!computed] / sqrt(2)))
  scaled.mean[!computed] = sqrt(2/pi) * expm1(-diff) / denom
  
  #double check that things are within bounds
  scaled.mean[scaled.mean > beta] = beta[scaled.mean > beta]
  scaled.mean[scaled.mean < alpha] = alpha[scaled.mean < alpha]
  
  # Flip back.
  scaled.mean[flip] = -scaled.mean[flip]
  
  #transform back to nonstandard normal case
  res[!sd.zero] = mean + sd * scaled.mean
  
  return(res)
}

#' @title Expected Squared Value of Truncated Normal
#' 
#' @description Computes the expected squared values of truncated normal 
#'   distributions with parameters \code{a}, \code{b}, \code{mean}, and 
#'   \code{sd}. Arguments can be scalars, vectors, or matrices. Arguments of 
#'   shorter length will be recycled according to the usual recycling rules, 
#'   but \code{a} and \code{b} must have the same length. Missing values are 
#'   accepted for all arguments.
#'
#' @inheritParams my_etruncnorm
#' 
#' @param sd The standard deviation of the untruncated normal. Standard
#'   deviations of zero are interpreted as numerically (rather than exactly)
#'   zero, so that the square of the untruncated mean is returned if it lies 
#'   within \code{[a, b]} and the square of the nearer of \code{a} and 
#'   \code{b} is returned otherwise.
#' 
#' @return The expected squared values of truncated normal
#'   distributions with parameters \code{a}, \code{b}, \code{mean}, and
#'   \code{sd}. If any of the arguments is a matrix, then a matrix will
#'   be returned.
#'  
#' @seealso \code{\link{my_etruncnorm}}, \code{\link{my_vtruncnorm}}
#'     
#' @export
#'
library(RcppFaddeeva)
my_e2truncnorm = function(a, b, mean = 0, sd = 1) {
  # based off of https://github.com/cossio/TruncatedNormal.jl
  do_truncnorm_argchecks(a, b)

  # initialize array to store result
  if (is.matrix(a)){
    res = array(dim = dim(0 * (a + b + mean + sd)))
  }
  else{
    res = rep(NA, length.out = length(0 * (a + b + mean + sd)))
  }

  #make sure input sizes all match 
  a = rep(a, length.out = length(res))
  b = rep(b, length.out = length(res))
  mean = rep(mean, length.out = length(res))
  sd = rep(sd, length.out = length(res))

  # Handle zero sd point masses
  sd.zero = (sd == 0)
  # if mean ≥ b, 2nd moment is b^2
  res[sd.zero & b <= mean] = b[sd.zero & b <= mean]^2
  # if mean ≤ a, 2nd moment is a^2
  res[sd.zero & a >= mean] = a[sd.zero & a >= mean]^2
  # if mean ∈ (a,b), 2nd moment is mean^2
  res[sd.zero & a < mean & mean < b] = mean[sd.zero & a < mean & mean < b]^2
  
  # Handle NAN inputs
  isna = is.na(a) | is.na(b) | is.na(mean) | is.na(sd)
  res[isna] = NA
  
  # Focus in on where sd is nonzero and nothing is nan
  a = a[!sd.zero & !isna]
  b = b[!sd.zero & !isna]
  mean = mean[!sd.zero & !isna]
  sd = sd[!sd.zero & !isna]

  # Rescale to standard normal distributions if sd is nonzero
  alpha = (a - mean) / sd
  beta = (b - mean) / sd
  scaled.mean = my_etruncnorm(alpha, beta)
  # initialize array for scaled 2nd moments
  scaled.2mom = rep(NA, length.out = length(alpha))

  # point mass bc endpoints are equal
  # 2nd moment is point-value squared
  endpoints_equal = (alpha == beta)
  scaled.2mom[endpoints_equal] = alpha[endpoints_equal] ^ 2
  # keep track of which spots in scaled.2mom are already computed
  computed = endpoints_equal

  # force to satisfy β ≥ 0 and |α| ≤ |β|
  # so either α ≤ 0 ≤ β or 0 < α ≤ β
  flip = !computed & abs(alpha) > abs(beta)
  flip[is.na(flip)] = FALSE
  orig.alpha = alpha
  alpha[flip] = -beta[flip]
  beta[flip] = -orig.alpha[flip]

  # both endpoints infinite/untruncated normal distribution
  # 2nd moment is 1
  both_inf = !computed & is.infinite(alpha) & is.infinite(beta)
  scaled.2mom[both_inf] = 1
  computed = computed | both_inf

  # truncated to [α,∞) 
  # 2nd moment simplifies to 1 + αϕ(α)/(1 - Φ(α))
  beta_inf = !computed & is.infinite(beta)
  scaled.2mom[beta_inf] = 1 + sqrt(2 / pi) * alpha[beta_inf] / Re(erfcx(alpha[beta_inf] / sqrt(2)))
  computed = computed | beta_inf

  # a ≤ 0 ≤ b
  #catestrophic cancellation is less of an issue
  alpha_negative = !computed & alpha <= 0
  ea = sqrt(pi/2) * Re(erf(alpha[alpha_negative] / sqrt(2)))
  eb = sqrt(pi/2) * Re(erf(beta[alpha_negative] / sqrt(2)))
  fa = ea - alpha[alpha_negative] * exp(-alpha[alpha_negative]^2 / 2)
  fb = eb - beta[alpha_negative] * exp(-beta[alpha_negative]^2 / 2)
  scaled.2mom[alpha_negative] = (fb - fa) / (eb - ea)
  computed = computed | alpha_negative

  # 0 < a ≤ b
  #strategically avoid catestrophic cancellation as much as possible
  exdiff = exp((alpha[!computed] - beta[!computed])*(alpha[!computed] + beta[!computed])/2)
  ea = sqrt(pi/2) * Re(erfcx(alpha[!computed] / sqrt(2)))
  eb = sqrt(pi/2) * Re(erfcx(beta[!computed] / sqrt(2)))
  fa = ea + alpha[!computed]
  fb = eb + beta[!computed]
  scaled.2mom[!computed] = (fa - fb * exdiff) / (ea - eb * exdiff)

  # transform results back to nonstandard normal case
  # μ^2 + σ^2 E(Z^2) + 2 μ σ E(Z)
  # possible catestropic cancellation because μ^2 + σ^2 E(Z^2) ≧ 0 while 2μσE(Z) has unknown sign
  # If |μ| < σ , compute as μ^2 + σ(σ E(Z^2) + 2 μ E(Z))
  mean_smaller = abs(mean) < sd
  res[!sd.zero & mean_smaller] = mean^2 + sd*(sd * scaled.2mom + 2 * mean * scaled.mean)
  # If  σ ≦ |μ| , compute as μ(μ +  2 σ E(Z)) + σ^2 E(Z^2)
  res[!sd.zero & !mean_smaller] = mean*(mean + 2 * sd * scaled.mean) + sd^2 * scaled.2mom

  # Check that the results make sense -- should never be negative
  stopifnot(all(res >= 0 | is.na(res)))  

  return(res)
}

#' @title Variance of Truncated Normal
#' @description Computes the variance of truncated normal distributions with
#'   parameters \code{a}, \code{b}, \code{mean}, and \code{sd}. Arguments can 
#'   be scalars, vectors, or matrices. Arguments of shorter length will be 
#'   recycled according to the usual recycling rules, but \code{a} and \code{b}
#'   must have the same length. Missing values are accepted for all arguments.
#'
#' @inheritParams my_etruncnorm
#' @param sd The standard deviation of the untruncated normal.
#' 
#' @return The variance of truncated normal distributions with parameters 
#'   \code{a}, \code{b}, \code{mean}, and \code{sd}. If any of the arguments 
#'   is a matrix, then a matrix will be returned.
#'   
#' @seealso \code{\link{my_etruncnorm}}, \code{\link{my_e2truncnorm}}
#'
#' @export
my_vtruncnorm = function(a, b, mean = 0, sd = 1) {
  # based off of https://github.com/cossio/TruncatedNormal.jl
  do_truncnorm_argchecks(a, b)
  
  #solved scaled problem
  alpha = (a - mean) / sd
  beta = (b - mean) / sd
  
  m1 = my_etruncnorm(alpha, beta)
  m2 = sqrt(my_e2truncnorm(alpha, beta))
  scaled.res = (m2 - m1) * (m1 + m2)
  
  # Handle endpoints equal
  scaled.res[(alpha == beta) & sd != 0] = alpha[(alpha == beta) & sd != 0]
  
  #transform back to unscaled
  res = sd^2 * scaled.res
  
  # Handle zero sds.
  res[sd == 0] = 0
  
  return(res)
}

do_truncnorm_argchecks = function(a, b) {
  if (!(length(a) == length(b)))
    stop("truncnorm functions require that a and b have the same length.")
  if (any(b < a, na.rm = TRUE))
    stop("truncnorm functions require that a <= b.")
}
