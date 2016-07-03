#sets optimization method
#deals with old flags VB and cxx
#also checks if necessary tools installed for optmethod specified
set_optmethod = function(optmethod,VB,cxx){
  # if user tries to set both optmethod and VB/cxx that's an error
  if ( !is.null(optmethod) && (!is.null(VB) || !is.null(cxx)) )
    stop("VB and cxx options are deprecated and incompatible with optmethod; use optmethod instead")
  
  # Fallbacks for optmethod
  # By default it will be "mixIP", if REBayes not present then fallback to EM
  if (!requireNamespace("REBayes", quietly = TRUE)) {  # check whether REBayes package is present
    # If REBayes package missing
    message("Due to absence of package REBayes, switching to EM algorithm")
    if (requireNamespace("Rcpp")) {
      optmethod = "cxxMixSquarem"
    } else {
      optmethod = "mixEM"  # fallback if neither Rcpp or REBayes are installed
      message("Using vanilla EM; for faster performance install REBayes (preferred) or Rcpp")
    }
  }
  
  # Check if VB and cxx are set to logical; for backwards compatibility
  if (!is.null(VB)) {
    warning("VB option is deprecated, use optmethod instead")
    if (VB == TRUE) optmethod = "mixVBEM"
  }
  
  if (!is.null(cxx)) {
    warning("cxx option is deprecated, use optmethod instead")
    if (cxx == TRUE)  optmethod = "cxxMixSquarem"
    if (cxx == FALSE) optmethod = "mixEM"
  }
  
  if (optmethod == "mixIP") assertthat::assert_that(requireNamespace("REBayes", quietly = TRUE))
  if (optmethod == "cxxMixSquarem") assertthat::assert_that(requireNamespace("Rcpp", quietly = TRUE))
  return(optmethod)
}


# handles defaults for optimization routines
set_control = function(control,nobs){
  # Handling control variables
  control.default = list(K = 1, method = 3, square = TRUE,
                         step.min0 = 1, step.max0 = 1, mstep = 4,
                         kr = 1, objfn.inc = 1, tol = 1.e-07, maxiter = 5000,
                         trace = FALSE)
  if (nobs > 50000) control.default$trace = TRUE
  namc = names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  control = modifyList(control.default, control)
  if (control$maxiter == 0)
    stop("option control$maxiter=0 deprecated; used fixg=TRUE instead")
  return(control)
}

check_args = function(mixcompdist,df,prior,optmethod,gridmult,sebetahat,betahat){
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



#' Takes raw data and sets up data object for use by ash
#' 
#' @details The data object stores both the data, and details of the model to be used for the data.
#' For example, in the generalized version of ash the cdf and pdf of the likelihood are
#' stored here.
#' 
#' @param betahat vector of betahats
#' @param sebetahat vector of standard errors
#' @param df degree of freedom to assume (for t likelihood)
#' @param alpha specifies value of alpha to use (model is for betahat/sebetahat^alpha | sebetahat)
#' 
#' @return data object (list) 
#' @export
set_data = function(betahat, sebetahat, df, alpha=0){
  # Dealing with precise input of betahat, currently we exclude them from the EM algorithm
  if(length(sebetahat)==1L){sebetahat = rep(sebetahat, length(betahat))}
  exclude = (sebetahat==0 | sebetahat == Inf | is.na(betahat) | is.na(sebetahat))
  
  data=list()
  data$x = betahat[!exclude]/(sebetahat[!exclude]^alpha)
  data$s = sebetahat[!exclude]^(1-alpha)
  data$v = df
  data$exclude = exclude
  data$alpha=alpha
  data$s_orig = sebetahat[!exclude]
  
  if(is.null(df)){
    data$cdfFUN = pnorm
    data$pdfFUN = dnorm
    data$etruncFUN = my_etruncnorm
    data$e2truncFUN = my_e2truncnorm
    data$FUNargs = list(mean = 0, sd = 1)
  } else {
    data$cdfFUN = ptgen
    data$pdfFUN = dtgen
    data$etruncFUN = my_etrunct
    data$e2truncFUN = my_e2trunct
    data$FUNargs = list(mean = 0, sd = 1, df=df)
  }
  
  return(data)
}