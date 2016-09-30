# This file contains functions written to compute 
# truncated expectations when no etruncFUN is provided


#' @title gen_etruncFUN
#' @description Produce function to compute expectation of truncated 
#'    error distribution from log cdf and log pdf (using numerical integration)
#'
#' @param lcdfFUN the log cdfFUN of the error distribution
#' @param lpdfFUN the log pdfFUN of the error distribution
gen_etruncFUN = function(lcdfFUN,lpdfFUN){
  return(function(a,b){
    tmp=mapply(gen_etruncFUN_single(lcdfFUN,lpdfFUN),a,b)
    dim(tmp) = dim(a)
    return(tmp)
  })
}


# compute expectation of truncated error distribution
# for scalars a and b
gen_etruncFUN_single = function(lcdfFUN,lpdfFUN){
  return(function(a,b){
    if(a == b){
      return(a)
    }else{
      denom = exp(lcdfFUN(b))-exp(lcdfFUN(a))
      if(denom!=0){ 
    
      # numerical integration
        xpdf = function(x){
          x*exp(lpdfFUN(x))
        }
        tmp = try(stats::integrate(xpdf,a,b)$value,silent=TRUE)
        if (class(tmp)!="try-error") 
          return(tmp/denom)
      }
    }
    return(NA)
  })
}

