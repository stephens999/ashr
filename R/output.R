# this file contains functions to compute the columns of
# the result data.frame from ash
# all such functions must take parameters (g,data)

# Posterior Mean
calc_pm = function(g,data){
  exclude = get_exclusions(data)
  PosteriorMean = rep(0,length = n_obs(data))
  PosteriorMean[!exclude] = postmean(g,data)[!exclude]
  PosteriorMean[exclude] = calc_mixmean(g)
  PosteriorMean = PosteriorMean * (data$s_orig^data$alpha)
  return(PosteriorMean)
}

# Posterior Standard Deviation
calc_psd = function(g,data){
  exclude  = get_exclusions(data)
  PosteriorSD = rep(0,length = n_obs(data))
  PosteriorSD[!exclude] = postsd(g,data)[!exclude]
  PosteriorSD[exclude] = calc_mixsd(g)
  PosteriorSD= PosteriorSD * (data$s_orig^data$alpha)
  return(PosteriorSD)
}

# Local FDR
calc_lfdr = function(g,data){
  exclude  = get_exclusions(data)
  ZeroProb = rep(0,length = n_obs(data))
  ZeroProb[!exclude] = colSums(comp_postprob(g,data)[pm_on_zero(g),,drop = FALSE])[!exclude]
  ZeroProb[exclude] = sum(mixprop(g)[pm_on_zero(g)])
  return(ZeroProb)
}

#negative probability
calc_np = function(g,data){
  exclude  = get_exclusions(data)
  NegativeProb = rep(0,length = n_obs(data))
  NegativeProb[!exclude] = cdf_post(g, 0, data)[!exclude] - calc_lfdr(g,data)[!exclude]
  NegativeProb[exclude] = mixcdf(g,0)
  return(NegativeProb)
}

#positive probability
calc_pp = function(g,data){
  pp=(1-calc_np(g,data)-calc_lfdr(g,data))
  ifelse(pp<0,0,pp) #deal with numerical issues that lead to numbers <0
}

# local False Sign Rate
calc_lfsr = function(g,data){
  compute_lfsr(calc_np(g,data),calc_lfdr(g,data))
}

calc_svalue = function(g,data){
  return(qval.from.lfdr(calc_lfsr(g,data)))
}

calc_qvalue = function(g,data){
  return(qval.from.lfdr(calc_lfdr(g,data)))
}

# Data for flashr package
calc_flash_data = function(g, data, penloglik) {
  kk = ncomp(g)
  n = n_obs(data)
  exclude = get_exclusions(data)
  comp_postprob = matrix(0, nrow = kk, ncol = n)
  comp_postmean = matrix(0, nrow = kk, ncol = n)
  comp_postmean2 =  matrix(0, nrow = kk, ncol = n)
  
  comp_postprob[, !exclude] = comp_postprob(g, data)[, !exclude]
  comp_postmean[, !exclude] = comp_postmean(g, data)[, !exclude]
  comp_postmean2[, !exclude] = comp_postmean2(g, data)[, !exclude]
  
  # For missing observations, use the prior instead of the posterior.
  comp_postprob[, exclude] = mixprop(g)
  comp_postmean[, exclude] = comp_mean(g)
  comp_postmean2[, exclude] = comp_mean2(g)
  
  postmean = colSums(comp_postprob * comp_postmean)
  postmean2 = colSums(comp_postprob * comp_postmean2)
  # Avoid potential negatives due to numeric rounding errors.
  postmean2[postmean2 < 0] = 0 
  
  return(list(fitted_g = g,
              postmean = postmean,
              postmean2 = postmean2,
              penloglik = penloglik))
}

# if outputlevel an integer, produce a vector of character strings naming output to be produced
# (otherwise return outputlevel)
set_output=function(outputlevel){
  if(!is.numeric(outputlevel)){return(outputlevel)} #allows that user might specify directly
  else{
    output = c("fitted_g")
    if(outputlevel>0){
      output = c(output,"loglik","logLR","PosteriorMean","PosteriorSD") 
    }
    if(outputlevel>1){output=c("data","NegativeProb","PositiveProb","lfsr","svalue","lfdr","qvalue",output)}
    if(outputlevel>2){output=c(output,"fit_details")}
    
    # These are special flags for output used by flashr.
    if(outputlevel==4){output=c("fitted_g","PosteriorMean", "PosteriorSD","flash_data")}
    if(outputlevel==5){output=c("flash_data")}
    
    return(output)
  }
}

# returns a named list of output functions, whose names are given in outputnames
# eg set_resfns(c("lfsr","lfdr")) would extract results functions for lfsr and lfdr
set_resfns = function(outputnames){
  result_fns = list(NegativeProb= calc_np, PositiveProb= calc_pp, lfsr = calc_lfsr, 
                  svalue = calc_svalue, lfdr = calc_lfdr, qvalue = calc_qvalue,PosteriorMean = calc_pm, PosteriorSD = calc_psd)

  return(result_fns[intersect(names(result_fns), outputnames)]) #extract the results functions specified in output 
}
