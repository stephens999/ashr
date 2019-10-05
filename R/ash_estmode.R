# a wrapper function that estimates the mode, using optimize
# called by ash if mode="estimate"
ash.estmode = function(betahat, modemin, modemax, ...){
  opt.fn = function(c) {
    return(-ash(betahat = betahat, mode = c, outputlevel = "loglik", ...)$loglik)
  }
  opt = stats::optimize(opt.fn, 
                        interval = c(modemin, modemax),
                        tol = abs(modemax - modemin) * .Machine$double.eps^0.25)
  return(opt$minimum)
}
