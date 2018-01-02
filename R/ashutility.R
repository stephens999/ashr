
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
#' @method summary ash
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
#' @method print ash
#' 
print.ash =function(x,...){
  print(summary(x,...))
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
#' @method plot ash
#' @export
#'
plot.ash = function(x,...,xmin,xmax){
  if(missing(xmin)){xmin=min(x$data$betahat)}
  if(missing(xmax)){xmax=max(x$data$betahat)}
  xgrid = seq(xmin,xmax,length=1000)
  y = cdf.ash(x,xgrid)
  graphics::plot(y,type="l",...)
}

#' @title Diagnostic plots for ash object
#'
#' @description Generate several plots to diagnose the fitness of ASH on the data
#'
#' @param x the fitted ash object
#' @param plot.it logical. whether to plot the diagnostic result
#' @param sebetahat.tol tolerance to test the equality of betahat
#' @param plot.hist logical. whether to plot the histogram of betahat when sebetahat is not constant
#' @param xmin,xmax range of the histogram of betahat to be plotted
#' @param breaks histograms parameter (see \code{\link[graphics]{hist}})
#' @param alpha error level for the de-trended diagnostic plot
#' @param pch,cex plot parameters for dots
#' 
#' @details None.
#'
#' @export
#'
#' @importFrom graphics lines
#' @importFrom graphics legend
#' @importFrom graphics abline
#' @importFrom stats punif
#' @importFrom stats qbeta
#' 
plot_diagnostic = function (x, plot.it = TRUE, 
                            sebetahat.tol = 1e-3,
                            plot.hist,
                            xmin, xmax, breaks = "Sturges",
                            alpha = 0.01,
                            pch = 19, cex = 0.25
                            ) {
  cdfhat = cdf_conv(x$fitted_g, x$data)
  na.ind = is.na(cdfhat)
  n = length(cdfhat[!na.ind])
  if (n == 0) (stop("The data have only NAs."))
  p.ks.unif = round(stats::ks.test(cdfhat, punif)$p.val, 3)
  upper = qbeta(1 - alpha / 2, 1:n, n + 1 - (1:n)) - (1:n) / (n + 1)
  lower = qbeta(alpha / 2, 1:n, n + 1 - (1:n)) - (1:n) / (n + 1)
  diff = sort(cdfhat[!na.ind]) - (1 : n) / (n + 1)
  if (plot.it) {
    sebetahat <- x$data$s
    sebetahat.same <- abs(max(sebetahat, na.rm = TRUE) - min(sebetahat, na.rm = TRUE)) / mean(sebetahat, na.rm = TRUE) <= sebetahat.tol
    if (missing(plot.hist)) {
      plot.hist = sebetahat.same
    }
    if (plot.hist) {
      betahat <- x$data$x
      if (missing(xmin)) {xmin = min(betahat, na.rm = TRUE)}
      if (missing(xmax)) {xmax = max(betahat, na.rm = TRUE)}
      xgrid.length = 1000
      xgrid = seq(xmin - 1, xmax + 1, length = xgrid.length)
      plot.data <- x$data
      if (sebetahat.same) {
        plot.data$x = xgrid
        plot.data$s = rep(mean(sebetahat, na.rm = TRUE), xgrid.length)
        fhat = dens_conv(x$fitted_g, plot.data)
      } else {
        fhat = c()
        for (i in 1 : xgrid.length) {
          plot.data$x = rep(xgrid[i], n)
          plot.data$s = sebetahat[!na.ind]
          fhat[i] = mean(dens_conv(x$fitted_g, plot.data))
        }
      }
      hist.betahat = graphics::hist(betahat[!na.ind], breaks = breaks, plot = FALSE)
      graphics::hist(betahat[!na.ind], probability = TRUE, breaks = breaks,
                     ylim = c(0, max(c(fhat, hist.betahat$density))),
                     xlab = expression(hat(beta)),
                     main = expression(paste("Histogram of ", hat(beta)))
                     )
      lines(xgrid, fhat, col = "blue")
      legend("topleft", lty = 1, col = "blue", "ASH")
      cat ("Press [enter] to see next plot")
      line <- readline()
    }
    graphics::plot((1 : n) / (n + 1), sort(cdfhat[!na.ind]),
                   xlim = c(0, 1), ylim = c(0, 1),
                   xlab = "Theoretical Uniform Quantile", 
                   ylab = "Estimated Predictive Quantile",
                   main = "Diagnostic Plot for ASH",
                   pch = pch, cex = cex
                   )
    abline(0, 1, lty = 2, col = "red")
    cat ("Press [enter] to see next plot")
    line <- readline()
    graphics::plot(diff, cex = cex, pch = pch, 
         ylim = range(diff, upper, lower, 0),
         xlab = "Index k",
         ylab = expression(Q[(k)] - E(Q[(k)])), 
         main = c("De-trended Diagnostic Plot for ASH", 
                  paste("K-S Uniformity Test p Value:", p.ks.unif))
         )
    abline(h = 0, col = "red", lty = 2)
    lines(upper, col = "red")
    lines(lower, col = "red")
    cat ("Press [enter] to see next plot")
    line <- readline()
    graphics::hist(cdfhat[!na.ind], probability = TRUE, breaks = breaks,
                   xlab = "Estimated Predictive Quantile",
                   main = c("Histogram of Estimated Predictive Quantile",
                            paste("K-S Uniformity Test p Value:", p.ks.unif))
                   )
    abline(h = 1, lty = 2, col = "red")
  }
  invisible(cdfhat)
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
  
  # compute log(dens_conv(g,data))
  log_comp_dens = log_comp_dens_conv(g, data)
  offset = apply(t(log_comp_dens),1,max)
  log_comp_dens = t(t(log_comp_dens)-offset) # avoid numeric issues by subtracting max of each row
  log_dens = log(colSums(g$pi * exp(log_comp_dens)))+offset # add offset back
  
  #return(log(dens_conv(g,data))- data$alpha*(log(data$s_orig)))
  return(log_dens - data$alpha*(log(data$s_orig)))
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

#' @title Sample from posterior 
#'
#' @description Returns random samples from the posterior distribution for each
#'     observation in an ash object. A matrix is returned, with columns corresponding
#'     to observations and rows corresponding to samples.
#'
#' @param a the fitted ash object
#' @param nsamp number of samples to return (for each observation)
#' @examples 
#' beta = rnorm(100,0,1)
#' betahat= beta+rnorm(100,0,1)
#' ash.beta = ash(betahat,1,mixcompdist="normal")
#' post.beta = get_post_sample(ash.beta,1000)
#' @export
get_post_sample = function(a,nsamp){
  return(post_sample(a$fitted_g,a$data,nsamp))
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





