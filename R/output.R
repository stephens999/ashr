# this file contains functions to compute the columns of
# the result data.frame from ash
# all such functions must take parameters (g,data)

# Posterior Mean
calc_pm = function(g,data){
  exclude = data$exclude
  PosteriorMean = rep(0,length = n_obs(data))
  PosteriorMean[!exclude] = postmean(g,data)
  PosteriorMean[exclude] = calc_mixmean(g)
  PosteriorMean = PosteriorMean * (data$sebetahat^data$alpha)
  return(PosteriorMean)
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

# Local FDR
calc_lfdr = function(g,data){
  exclude  = data$exclude
  ZeroProb = rep(0,length = n_obs(data))
  ZeroProb[!exclude] = colSums(comp_postprob(g,data)[pm_on_zero(g),,drop = FALSE])
  ZeroProb[exclude] = sum(mixprop(g)[pm_on_zero(g)])
  return(ZeroProb)
}

#negative probability
calc_np = function(g,data){
  exclude  = data$exclude
  NegativeProb = rep(0,length = n_obs(data))
  NegativeProb[!exclude] = cdf_post(g, 0, data) - calc_lfdr(g,data)[!exclude]
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


# Data for flash package
calc_flash_data = function(g,data){
  kk = ncomp(g)
  n = n_obs(data)
  exclude = data$exclude
  comp_postprob = matrix(0,nrow = kk, ncol = n)
  comp_postmean = matrix(0,nrow = kk, ncol = n)
  comp_postmean2 =  matrix(0,nrow = kk, ncol = n)
  
  comp_postprob[,!exclude] = comp_postprob(g,data)
  comp_postmean[,!exclude] = comp_postmean(g,data)
  comp_postmean2[,!exclude] = comp_postmean2(g,data)
  
  #FOR MISSING OBSERVATIONS, USE THE PRIOR INSTEAD OF THE POSTERIOR
  comp_postprob[,exclude] = mixprop(g)
  comp_postmean[,exclude] = comp_mean(g)
  comp_postmean2[,exclude] = comp_mean2(g)
  return(list(comp_postprob = comp_postprob,comp_postmean = comp_postmean,comp_postmean2 = comp_postmean2))
}


# If outputlevel is a list, then just returns it
# if outputlevel an integer, there are different combinations of
# default output provided
set_output=function(outputlevel){
  if(is.list(outputlevel)){output=outputlevel} 
  else {
    output = list(fitted.g=TRUE, call=TRUE)
    if(outputlevel>0){output$loglik=TRUE; output$logLR=TRUE
    output$resfns = list(PosteriorMean = calc_pm, PosteriorSD = calc_psd)
    }
    if(outputlevel>1){output$data=TRUE
      output$resfns = c(NegativeProb = calc_np, PositiveProb= calc_pp, lfsr = calc_lfsr, svalue = calc_svalue, 
                      lfdr = calc_lfdr, qvalue = calc_qvalue, output$resfns)
    }
    if(outputlevel>2){output$fit_details=TRUE}
    if(outputlevel>3){output$flash.data = TRUE}
  }
  return(output)
}
