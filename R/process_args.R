# avoid "no visible binding for global variable" note in CRAN check
# These variables are actually defined in process_args
if(getRversion() >= "2.15.1") utils::globalVariables(c("VB","cxx","method","model","mixcompdist","gridmult","control"))

#' Process input arguments for ash.workhorse
#'
#' @param oldargs captured argument list
#' @return list containing the processed arguments
#' @importFrom utils modifyList
#' @export process_args
process_args = function (oldargs) {
  # Assign each captured argument in the list to a variable
  for (i in 1L:length(oldargs)) assign(names(oldargs)[i], oldargs[[i]])
  
  # Start processing arguments
  
  if (length(sebetahat) == 1L) sebetahat = rep(sebetahat, length(betahat))
  if (length(sebetahat) != length(betahat))
    stop("Error: sebetahat must have length 1, or same length as betahat")
  
  # Set optimization method (optmethod)
  
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
  
  # method provides a convenient interface to set a particular combinations of parameters for prior an
  # If method is supplied, use it to set up specific values for these parameters; provide warning if values
  # are also specified by user
  # If method is supplied use the user-supplied values (or defaults if user does not specify them)
  if (method == "shrink") {
    
    # Almost equivalent to is.missing(prior)
    if (identical(sort(prior), sort(c("nullbiased","uniform","unit")))) {
      prior = "uniform"
    } else {
      warning("Specification of prior overrides default for method shrink")
    }
    
    if (is.null(pointmass)) {
      pointmass = FALSE
    } else if (pointmass != FALSE) {
      warning("Specification of pointmass overrides default for method shrink")
    }
    
  }
  
  if (method == "fdr") {
    
    # Almost equivalent to is.missing(prior)
    if (identical(sort(prior), sort(c("nullbiased","uniform","unit")))) {
      prior = "nullbiased"
    } else {
      warning("Specification of prior overrides default for method fdr")
    }
    
    if (is.null(pointmass)) {
      pointmass = TRUE
    } else if (pointmass != TRUE) {
      warning("Specification of pointmass overrides default for method fdr")
    }
    
  }
  
  # Dealing with precise input of betahat, currently we exclude them from the EM algorithm
  betahat.input = betahat
  sebetahat.input = sebetahat
  excludeindex = c(1:length(sebetahat.input))[sebetahat.input==0]
  if(length(excludeindex) == 0) excludeindex = NULL
  betahat = betahat.input[sebetahat.input != 0]
  sebetahat = sebetahat.input[sebetahat.input != 0]
  
  # Set observations with infinite standard errors to missing
  # later these missing observations will be ignored in EM, and posterior will be same as prior.
  sebetahat[sebetahat == Inf] = NA
  betahat[sebetahat == Inf] = NA
  
  if (model == "ET") {  # for ET model, standardize
    betahat = betahat/sebetahat
    sebetahat.orig = sebetahat  # store so that can be reinstated later
    sebetahat = rep(1, length(betahat))
  }
  
  if (mixcompdist == "normal" & !is.null(df))
    stop("Error: Normal mixture for student-t likelihood is not yet implemented")
  
  if (identical(prior, "unit") & optmethod != "mixVBEM")
    stop("Error: unit prior only valid for mixVBEM")
  
  if (mixcompdist == "halfuniform" & !identical(prior, "nullbiased"))
    warning("Use of halfuniform without nullbiased prior can lead to misleading local false sign rates, and so is not recommended")
  
  if (gridmult <= 1) stop("gridmult must be > 1")
  
  completeobs = (!is.na(betahat) & !is.na(sebetahat))
  n = sum(completeobs)
  
  # Handling control variables
  control.default = list(K = 1, method = 3, square = TRUE,
                         step.min0 = 1, step.max0 = 1, mstep = 4,
                         kr = 1, objfn.inc = 1, tol = 1.e-07, maxiter = 5000,
                         trace = FALSE)
  if (n > 50000) control.default$trace = TRUE
  namc = names(control)
  if (!all(namc %in% names(control.default)))
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  controlinput = modifyList(control.default, control)
  if (controlinput$maxiter == 0)
    stop("option control$maxiter=0 deprecated; used fixg=TRUE instead")
  
  if (n == 0) stop("Error: all input values are missing")
  
  # Collect everything into a new list
  newargs_names = setdiff(ls(), c("oldargs", "i", "call"))
  newargs = list()
  # assigning NULL to list component will remove that component; so use lapply
  safe_assign = function(x) if (is.null(get(x))) return(NULL) else return(get(x))
  newargs = lapply(newargs_names, safe_assign)
  names(newargs) = newargs_names
  
  return(newargs)
  
}
