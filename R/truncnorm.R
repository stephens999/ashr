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
#' 
my_etruncnorm = function(a, b, mean = 0, sd = 1) {
  do_truncnorm_argchecks(a, b)
  
  # The case where some sds are zero is handled last. In the meantime, assume
  #   that sd > 0.
  alpha = (a - mean) / sd
  beta = (b - mean) / sd
  
  # Flip alpha and beta when: 1. Both are positive (since computations are 
  #   unstable when both values of pnorm are close to 1); 2. dnorm(alpha) is 
  #   greater than dnorm(beta) (since subtraction is done on the log scale).
  flip = (alpha > 0 & beta > 0) | (beta > abs(alpha))
  flip[is.na(flip)] = FALSE
  orig.alpha = alpha
  alpha[flip] = -beta[flip]
  beta[flip] = -orig.alpha[flip]
  
  dnorm.diff = logscale_sub(dnorm(beta, log = TRUE), dnorm(alpha, log = TRUE))
  pnorm.diff = logscale_sub(pnorm(beta, log.p = TRUE), pnorm(alpha, log.p = TRUE))
  scaled.res = -exp(dnorm.diff - pnorm.diff)
  
  # Handle the division by zero that occurs when pnorm.diff = -Inf (that is, 
  #   when endpoints are approximately equal).
  endpts.equal = is.infinite(pnorm.diff)
  scaled.res[endpts.equal] = (alpha[endpts.equal] + beta[endpts.equal]) / 2
  
  # When alpha and beta are very large and both negative (due to the flipping
  #   logic, they cannot both be positive), computations can become unstable. 
  #   We find such cases by checking that the expectations make sense. When
  #   beta is negative, beta + 1 / beta is a lower bound for the expectation.
  #   Further, it is an increasingly good approximation as beta goes to -Inf 
  #   as long as alpha and beta aren't too close to one another. If they are,
  #   then their midpoint can be used as an alternate approximation (and lower 
  #   bound).
  lower.bd = pmax(beta + 1 / beta, (alpha + beta) / 2)
  bad.idx = (!is.na(beta) & beta < 0 
             & (scaled.res < lower.bd | scaled.res > beta))
  scaled.res[bad.idx] = lower.bd[bad.idx]
  
  # Flip back.
  scaled.res[flip] = -scaled.res[flip]
  
  res = mean + sd * scaled.res
  
  # Handle zero sds. Return the mean of the untruncated normal when it is 
  #   located inside of the interval [alpha, beta]. Otherwise, return the 
  #   endpoint that is closer to the mean.
  if (any(sd == 0)) {
    # For the subsetting to work correctly, arguments need to be recycled.
    a = rep(a, length.out = length(res))
    b = rep(b, length.out = length(res))
    mean = rep(mean, length.out = length(res))
    
    sd.zero = (sd == 0)
    res[sd.zero & b <= mean] = b[sd.zero & b <= mean]
    res[sd.zero & a >= mean] = a[sd.zero & a >= mean]
    res[sd.zero & a < mean & b > mean] = mean[sd.zero & a < mean & b > mean]
  }
  
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
    # based off of https://github.com/cossio/TruncatedNormal.jl/blob/3c16866c3afa3920e787513d492689e9e81192ca/src/tnmom2.jl
    do_truncnorm_argchecks(a, b)

    if ((mean == 0) && (sd ==1)){
        #standard normal distribution
        if (a == b){
            # point mass 
            # ⟹ 2nd moment is point-value squared
            return(a^2)
        } 
        else if (abs(a) > abs(b)) {
            # forces a and b to satisfy b ≥ 0, |a| ≤ |b|
            return(my_e2truncnorm(-b, -a))
        } 
        else if (is.infinite(a) && is.infinite(b)){
            # untruncated normal distribution
            # ⟹ 2nd moment is 1
            return(1.0)
        } 
        else if (is.infinite(b)) {
            # truncated to [a,∞) 
            # ⟹ 2nd moment simplifies to 1 + aϕ(a)/(1 - Φ(a))
            m2 = 1 + sqrt(2 / pi) * a / Re(erfcx(a / sqrt(2)))
            stopifnot(a^2 <= m2)
            return(m2)
        }
             
        #now everything is finite, b ≥ 0, and |a| ≤ |b|, so 
        # either a ≤ 0 ≤ b or 0 ≤ a ≤ b
        stopifnot(a < b, b < Inf, abs(a) <= abs(b))
        stopifnot(((a <= 0) && (0 <= b)) || ((0 <= a) && (a <= b)))

        if ((a <= 0) && (0 <= b)){ 
            # a ≤ 0 ≤ b
            #catestrophic cancellation is less of an issue
            ea = sqrt(pi/2) * Re(erf(a / sqrt(2)))
            eb = sqrt(pi/2) * Re(erf(b / sqrt(2)))
            fa = ea - a * exp(-a^2 / 2)
            fb = eb - b * exp(-b^2 / 2)
            m2 = (fb - fa) / (eb - ea)
            stopifnot(fb >= fa, eb >= ea)
            stopifnot(0 <= m2, m2 <= 1)
            return(m2)
        }  
        else { 
            # 0 ≤ a ≤ b
            #strategically avoid catestrophic cancellation as much as possible
            exdiff = exp((a - b)*(a + b)/2)
            ea = sqrt(pi/2) * Re(erfcx(a / sqrt(2)))
            eb = sqrt(pi/2) * Re(erfcx(b / sqrt(2)))
            fa = ea + a
            fb = eb + b
            m2 = (fa - fb * exdiff) / (ea - eb * exdiff)
            stopifnot(a^2 <= m2, m2 <= b^2)
            return(m2)
        }
    }
    else if (sd > 0){
        #nonstandard normal
        alpha = (a - mean) / sd
        beta = (b - mean) / sd
        # TODO potential for catestrophic cancellation here... is there a better way?
        return(mean^2 + sd^2 * my_e2truncnorm(alpha, beta) + 2 * mean * sd * my_etruncnorm(alpha, beta))
    } 
    else {
        #point mass
        # TODO potential error in original the julia code??
        # i wrote what I think it should be, but that doesn't agree with 
        # the original code
        if ((a <= mean) && (a <= b)) {
            # ⟹ if mean ∈ [a,b], 2nd moment is mean^2
            return(mean^2)
        } else if (mean < a){
            # ⟹ if mean < a, 2nd moment is a^2
            return(a^2)
        } else if (mean > b){
            # ⟹ if mean > a, 2nd moment is b^2
            return(b^2)
        }
    }
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
#' 
my_vtruncnorm = function(a, b, mean = 0, sd = 1) {
  do_truncnorm_argchecks(a, b)
  
  alpha = (a - mean) / sd
  beta = (b - mean) / sd
  
  scaled.res = my_e2truncnorm(alpha, beta) - my_etruncnorm(alpha, beta)^2
  
  # When alpha and beta are large and share the same sign, this computation 
  #   becomes unstable. A good approximation in this regime is 1 / beta^2
  #   (when alpha and beta are both negative). If my_e2truncnorm and
  #   my_etruncnorm are accurate to the eighth digit, then we can only trust 
  #   results up to the second digit if beta^2 and 1 / beta^2 differ by an
  #   order of magnitude no more than 6.
  smaller.endpt = pmin(abs(alpha), abs(beta))
  bad.idx = (is.finite(smaller.endpt) & smaller.endpt > 30)
  scaled.res[bad.idx] = pmin(1 / smaller.endpt[bad.idx]^2,
                             (beta[bad.idx] - alpha[bad.idx])^2 / 12)
  
  # Handle zero sds.
  scaled.res[is.nan(scaled.res)] = 0
  
  res = sd^2 * scaled.res
  
  return(res)
}

do_truncnorm_argchecks = function(a, b) {
  if (!(length(a) == length(b)))
    stop("truncnorm functions require that a and b have the same length.")
  if (any(b < a, na.rm = TRUE))
    stop("truncnorm functions require that a <= b.")
}

logscale_sub = function(logx, logy) {
  # In rare cases, logx can become numerically less than logy. When this
  #   occurs, logx is adjusted and a warning is issued.
  diff = logx - logy
  if (any(diff < 0, na.rm = TRUE)) {
    bad.idx = (diff < 0)
    bad.idx[is.na(bad.idx)] = FALSE
    logx[bad.idx] = logy[bad.idx]
    warning("logscale_sub encountered negative value(s) of logx - logy (min: ",
            formatC(min(diff[bad.idx]), format = "e", digits = 2), ")")
  }
  
  scale.by = logx
  scale.by[is.infinite(scale.by)] = 0
  return(log(exp(logx - scale.by) - exp(logy - scale.by)) + scale.by)
}
