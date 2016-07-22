# this file contains functions to compute outputs from ash
# The calc_xx functions just compute a result
# The output_xx functions return a named list, 
# with the names being used in the ash output

# Posterior Mean
calc_pm = function(g,data){
  exclude = data$exclude
  PosteriorMean = rep(0,length = n_obs(data))
  PosteriorMean[!exclude] = postmean(g,data)
  PosteriorMean[exclude] = calc_mixmean(g)
  PosteriorMean = PosteriorMean * (data$sebetahat^data$alpha)
  return(PosteriorMean)
}
output_pm = function(g,data){
  return(list(PosteriorMean = calc_pm(g,data)))  
}

# Posterior Standard Deviation
calc_psd = function(g,data){
  exclude  = data$exclude
  PosteriorSD = rep(0,length = n_obs(data))
  PosteriorSD[!exclude] = postsd(g,data)
  PosteriorSD[exclude] = calc_mixsd(g)
  PosteriorSD= PosteriorSD * (data$sebetahat^data$alpha)
  return(PosteriorSD)
}
output_psd = function(g,data){
  list(PosteriorSD = calc_psd(g,data))
}

# Local FDR
calc_lfdr = function(g,data){
  exclude  = data$exclude
  ZeroProb = rep(0,length = n_obs(data))
  ZeroProb[!exclude] = colSums(comppostprob(g,data)[comp_sd(g) ==0,,drop = FALSE])
  ZeroProb[exclude] = sum(mixprop(g)[comp_sd(g) == 0])
  return(ZeroProb)
}
output_lfdr = function(g,data){
  list(lfdr=calc_lfdr(g,data))
}

#negative probability
calc_np = function(g,data){
  exclude  = data$exclude
  NegativeProb = rep(0,length = n_obs(data))
  NegativeProb[!exclude] = cdf_post(g, 0, data) - calc_lfdr(g,data)[!exclude]
  NegativeProb[exclude] = mixcdf(g,0)
  return(NegativeProb)
}
output_np = function(g,data){
 list(NegativeProb = calc_np(g,data)) 
}

#positive probability
calc_pp = function(g,data){
  pp=(1-calc_np(g,data)-calc_lfdr(g,data))
  ifelse(pp<0,0,pp) #deal with numerical issues that lead to numbers <0
}
output_pp = function(g,data){
  list(PositiveProb = calc_pp(g,data))
}

# local False Sign Rate
calc_lfsr = function(g,data){
  compute_lfsr(calc_np(g,data),calc_lfdr(g,data))
}
output_lfsr = function(g,data){
  return(list(lfsr=calc_lfsr(g,data)))
}

output_svalue = function(g,data){
  return(list(svalue = qval.from.lfdr(calc_lfsr(g,data))))
}

output_qvalue = function(g,data){
  return(list(qvalue = qval.from.lfdr(calc_lfdr(g,data))))
}


# Data for flash package
calc_flash_data = function(g,data){
  kk = ncomp(g)
  n = n_obs(data)
  exclude = data$exclude
  comp_postprob = matrix(0,nrow = kk, ncol = n)
  comp_postmean = matrix(0,nrow = kk, ncol = n)
  comp_postmean2 =  matrix(0,nrow = kk, ncol = n)
  
  comp_postprob[,!exclude] = comppostprob(g,data)
  comp_postmean[,!exclude] = comp_postmean(g,data)
  comp_postmean2[,!exclude] = comp_postmean2(g,data)
  
  #FOR MISSING OBSERVATIONS, USE THE PRIOR INSTEAD OF THE POSTERIOR
  comp_postprob[,exclude] = mixprop(g)
  comp_postmean[,exclude] = comp_mean(g)
  comp_postmean2[,exclude] = comp_mean2(g)
  return(list(comp_postprob = comp_postprob,comp_postmean = comp_postmean,comp_postmean2 = comp_postmean2))
}
output_flash_data = function(g,data){
  return(list(flash.data=calc_flash_data(g,data)))
}


#outputs list for output
output_loglik = function(g,data){
  return(list(loglik = calc_loglik(g,data) ))
}

#outputs list for output
output_logLR = function(g,data){
  return(list(logLR = calc_logLR(g,data) ))
}