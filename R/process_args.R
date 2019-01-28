#sets optimization method
#also checks if necessary tools installed for optmethod specified
set_optmethod = function(optmethod){
  # Fallbacks for optmethod - checks required packages are installed
  if(optmethod == "mixIP"){
      if (!requireNamespace("REBayes", quietly = TRUE)) {
          
    # Check whether REBayes package is present.
    # If REBayes package missing.
      message("Due to absence of package REBayes, switching to EM algorithm")
      optmethod = "cxxMixSquarem"
    }
  }
  
  if (optmethod == "mixIP")
    assertthat::assert_that(requireNamespace("REBayes", quietly = TRUE))
  return(optmethod)
}

check_lik = function(lik, betahat, sebetahat, df, mixcompdist){
  if(is.null(lik$lcdfFUN)){stop("Likelihood must have lcdfFUN")}
  if(is.null(lik$lpdfFUN)){stop("Likelihood must have lpdfFUN")}
  
  if(!(lik$name %in% c("normal","t")) & !is.null(df)){
    warning("Input df is ignored for this likelihood function")
  }
  if(!(lik$name %in% c("normal","normalmix")) & mixcompdist == "normal"){
    stop("Error: Normal mixture for non-normal likelihood is not yet implemented")
  }
  if(lik$name %in% c("pois","binom")){
    if (!sum(betahat==0) | !sum(sebetahat==1)){
      stop("Error: betahat must be rep(0,n) and sebetaht must be 1 for 
           Poisson/Binomial likelihood")
    }
  }
  
}


check_args = function(mixcompdist,df,prior,optmethod,gridmult,sebetahat,betahat){
  if(!is.numeric(betahat)){
    stop("Error: betahat must be numeric")
  }
  if(!is.numeric(sebetahat)){
    stop("Error: sebetahat must be numeric")
  }
  
  if (mixcompdist == "normal" & !is.null(df))
    stop("Error: Normal mixture for student-t likelihood is not yet implemented")
  
  if (identical(prior, "unit") & optmethod != "mixVBEM")
    stop("Error: unit prior only valid for mixVBEM")
  
  if (mixcompdist == "halfuniform" & !identical(prior, "nullbiased"))
    warning("Use of halfuniform without nullbiased prior can lead to misleading local false sign rates, and so is not recommended")
  
  if (gridmult <= 1) stop("gridmult must be > 1")
  
  if ((length(sebetahat) != length(betahat)) & (length(sebetahat) != 1))
    stop("Error: sebetahat must have length 1, or same length as betahat")
  
  completeobs = (!is.na(betahat) & !is.na(sebetahat))
  if (sum(completeobs) == 0) stop("Error: all input values are missing")
}



