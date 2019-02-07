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
my_e2truncnorm = function(a, b, mean = 0, sd = 1) {
  do_truncnorm_argchecks(a, b)
  
  alpha = (a - mean) / sd
  beta = (b - mean) / sd
  
  # Flip alpha and beta when both are positive (as above, but the mean is
  #   also recycled and flipped so that we don't have to flip back).
  flip = (alpha > 0 & beta > 0)
  flip[is.na(flip)] = FALSE
  orig.alpha = alpha
  alpha[flip] = -beta[flip]
  beta[flip] = -orig.alpha[flip]
  if (any(mean != 0)) {
    mean = rep(mean, length.out = length(alpha))
    mean[flip] = -mean[flip]
  }
  
  pnorm.diff = logscale_sub(pnorm(beta, log.p = TRUE), pnorm(alpha, log.p = TRUE))
  alpha.frac = alpha * exp(dnorm(alpha, log = TRUE) - pnorm.diff)
  beta.frac = beta * exp(dnorm(beta, log = TRUE) - pnorm.diff)
  
  # Create a vector or matrix of 1's with NA's in the correct places.
  if (is.matrix(alpha))
    scaled.res = array(1, dim = dim(alpha))
  else
    scaled.res = rep(1, length.out = length(alpha))
  is.na(scaled.res) = is.na(flip)
  
  alpha.idx = is.finite(alpha)
  scaled.res[alpha.idx] = 1 + alpha.frac[alpha.idx]
  beta.idx = is.finite(beta)
  scaled.res[beta.idx] = scaled.res[beta.idx] - beta.frac[beta.idx]
  
  # Handle approximately equal endpoints.
  endpts.equal = is.infinite(pnorm.diff)
  scaled.res[endpts.equal] = (alpha[endpts.equal] + beta[endpts.equal])^2 / 4
  
  # Check that the results make sense. When beta is negative,
  #   beta^2 + 2 * (1 + 1 / beta^2) is an upper bound for the expected squared
  #   value, and it is typically a good approximation as beta goes to -Inf. 
  #   When the endpoints are very close to one another, the expected squared 
  #   value of the uniform distribution on [alpha, beta] is a better upper 
  #   bound (and approximation).
  upper.bd1 = beta^2 + 2 * (1 + 1 / beta^2)
  upper.bd2 = (alpha^2 + alpha * beta + beta^2) / 3
  upper.bd = pmin(upper.bd1, upper.bd2)
  bad.idx = (!is.na(beta) & beta < 0 
             & (scaled.res < beta^2 | scaled.res > upper.bd))
  scaled.res[bad.idx] = upper.bd[bad.idx]
  
  res = mean^2 + 2 * mean * sd * my_etruncnorm(alpha, beta) + sd^2 * scaled.res
  
  # Handle zero sds.
  if (any(sd == 0)) {
    a = rep(a, length.out = length(res))
    b = rep(b, length.out = length(res))
    mean = rep(mean, length.out = length(res))
    
    sd.zero = (sd == 0)
    res[sd.zero & b <= mean] = b[sd.zero & b <= mean]^2
    res[sd.zero & a >= mean] = a[sd.zero & a >= mean]^2
    res[sd.zero & a < mean & b > mean] = mean[sd.zero & a < mean & b > mean]^2
  }
  
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
