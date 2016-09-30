
################################## ASH UTILITY FUNCTIONS ############################

#' @title Summary method for ash object
#'
#' @description Print summary of fitted ash object
#'
#' @details \code{\link{summary}} prints the fitted mixture, the
#'     fitted log likelihood with 10 digits and a flag to indicate
#'     convergence
#' @param object the fitted ash object
#' @param ... not used, included for consistency as an S3
#'     generic/method.
#'
#' @export
#'
summary.ash=function(object,...){
  if(!is.null(object$fitted_g)){print(object$fitted_g)}
  print("Ash object has elements: ")
  print(paste(names(object),sep=","))
}


#' @title Print method for ash object
#'
#' @description Print the fitted distribution of beta values in the EB
#'     hierarchical model
#'
#' @details None
#' @param x the fitted ash object
#' @param ... not used, included for consistency as an S3
#'     generic/method.
#'
#' @export
#'
print.ash =function(x,...){
  summary(x)
}

#' @title Plot method for ash object
#'
#' @description Plot the cdf of the underlying fitted distribution
#'
#' @param x the fitted ash object
#' @param ... Arguments to be passed to methods,such as graphical parameters (see \code{\link[graphics]{plot}})
#' @param xmin xlim lower range, default is the lowest value of betahat
#' @param xmax xlim upper range, default is the highest value of betahat
#' @details None
#'
#' @export
#'
plot.ash = function(x,...,xmin,xmax){
  if(missing(xmin)){xmin=min(x$data$betahat)}
  if(missing(xmax)){xmax=max(x$data$betahat)}
  xgrid = seq(xmin,xmax,length=1000)
  y = cdf.ash(x,xgrid)
  graphics::plot(y,type="l",...)
}

#' @title Compute loglikelihood for data from ash fit
#'
#' @description Return the log-likelihood of the data for a given g() prior 
#'
#' @param g the fitted g, or an ash object containing g
#' @param data a data object, see set_data
#'
#' @export
calc_loglik = function(g,data){
  sum(calc_vloglik(g,data))
}

#' @title Compute loglikelihood for data under null that all beta are 0
#'
#' @description Return the log-likelihood of the data betahat, with
#'     standard errors betahatsd, under the null that beta==0
#'
#' @param data a data object; see set_data
#'
#' @export
calc_null_loglik = function(data){
  sum(calc_null_vloglik(data))
}


#' @title Compute loglikelihood ratio for data from ash fit
#'
#' @description Return the log-likelihood ratio of the data for a given g() prior 
#'
#' @inheritParams calc_loglik
#' @export
calc_logLR = function(g,data){
  return(calc_loglik(g,data) - calc_null_loglik(data))
}

#' @title Compute vector of loglikelihood for data from ash fit
#'
#' @description Return the vector of log-likelihoods of the data
#'     betahat, with standard errors betahatsd, for a given g() prior
#'     on beta, or an ash object containing that
#'
#' @inheritParams calc_loglik
#'
#' @export
#'
calc_vloglik = function(g,data){
  if(class(g)=="ash"){
    if(g$data$alpha != data$alpha){
      warning("Model (alpha value) used to fit ash does not match alpha in data! Probably you have made a mistake!")
    }
    if(class(g)=="ash"){g = g$fitted_g} #extract g object from ash object if ash object passed
  }
  return(log(dens_conv(g,data))- data$alpha*(log(data$s_orig)))
}


#' @title Compute vector of loglikelihood for data under null that all
#'     beta are 0
#'
#' @description Return the vector of log-likelihoods of the data points under the null
#'
#' @inheritParams calc_null_loglik
#'
#' @export
calc_null_vloglik = function(data){
    return(do.call(data$lik$lpdfFUN, list(x=data$x/data$s)) - log(data$s)
           -data$alpha*log(data$s_orig))
}

#' @title Compute vector of loglikelihood ratio for data from ash fit
#'
#' @description Return the vector of log-likelihood ratios of the data
#'     betahat, with standard errors betahatsd, for a given g() prior
#'     on beta, or an ash object containing that, vs the null that g()
#'     is point mass on 0
#'
#' @inheritParams calc_loglik
#'
#' @export
calc_vlogLR = function(g,data){
  return(calc_vloglik(g,data) - calc_null_vloglik(data))
}


#' @title Density method for ash object
#'
#' @description Return the density of the underlying fitted distribution
#'
#' @param a the fitted ash object
#' @param x the vector of locations at which density is to be computed
#'
#' @details None
#'
#' @export
#'
#'
get_density=function(a,x){
  list(x=x,y=dens(a$fitted_g,x))
}


#' @title cdf method for ash object
#'
#' @description Computed the cdf of the underlying fitted distribution
#'
#' @param a the fitted ash object
#' @param x the vector of locations at which cdf is to be computed
#' @param lower.tail (default=TRUE) whether to compute the lower or upper tail
#'
#' @details None
#'
#' @export
cdf.ash=function(a,x,lower.tail=TRUE){
  return(list(x=x,y=mixcdf(a$fitted_g,x,lower.tail)))
}


#Functions from MATLAB packages, used to measure performance and to show progress
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]
  assign(".tic", tic, envir=baseenv())
  invisible(tic)
}

toc <- function()
{
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  print(toc - tic)
  invisible(toc)
}





