
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




#' @title Multi-model Adaptive Shrinkage function
#'
#' @description This is a wrapper function that takes a grid value of
#'     \eqn{alpha} and then consider the model \eqn{betahat_j /
#'     s_j^{alpha} ~ g()},and eqn{beta_j / s_j^{alpha} ~
#'     N(0,(sebetahat^(1-alpha))^2) or student t distribution}. When
#'     \eqn{alpha=0} we have the EE model, when \eqn{alpha=1}, we have
#'     the ET model. \eqn{alpha} should be in between 0 and 1,
#'     inclusively. This wrapper function would select the best
#'     \eqn{alpha} and reports the ash item based on that \eqn{alpha}.
#'
#' @seealso \code{\link{ash}} the main function that this wrapper
#'     function is calling
#' @details All other inputs are exactly the same as the main function
#'     ash, and would pass to the main function to evaluate the
#'     likelihood.
#'
#' @param betahat  a p vector of estimates
#' @param sebetahat a p vector of corresponding standard errors

#' @param mixcompdist distribution of components in mixture. Default
#'     is "uniform".
#' @param df appropriate degrees of freedom for (t) distribution of
#'     betahat/sebetahat, default is NULL(Gaussian)
#' @param alpha Could be a vector of grid values in interval [0,1],
#'     that this wrapper would select based on likelihood
#'     principle. Could also be a positive integer greater or equal to
#'     2, then alpha number of grid values would be generated from
#'     [0,1], equally spaced. The default value is 2 that we compare
#'     the EE and ET model.
#' @param ncores Whether to use parallel computing, defaults to FALSE,
#'     user could specify number of cores they would like to
#'     use. Further, if user does not specify and
#'     length(betahat)>50000, then the function would perform parallel
#'     computation using number of CPU cores on the current host.
#' @param ... Further arguments to be passed to \code{\link{ash}}.
#'
#' @return ashm returns a list of objects
#' \item{beta.ash}{the best fitted ash object}
#' \item{loglikvector}{the vector of loglikelihood of various models}
#' \item{allash}{the fitted ash of various models}
#'

#' @export
#' @examples
#' beta = c(rep(0,100),rnorm(100))
#' sebetahat = abs(rnorm(200,0,1))
#' betahat = rnorm(200,beta,sebetahat)
#' beta.ashm = ashm(betahat, sebetahat,alpha=6)
#' #beta.ashm4 = ashm(betahat, sebetahat,alpha=6,ncores=4)
#' print(beta.ashm[[1]])  #best ash object
#' print(beta.ashm[[2]])  #corresponding model type
#' print(beta.ashm[[3]])  #log-likelihood for all models
#'
#'
ashm=function(betahat,sebetahat,
              mixcompdist = c("uniform","halfuniform","normal","+uniform","-uniform"),
			  df=NULL, alpha=2,ncores=FALSE,
              ...){
  if(length(alpha)==1){
    alpha=seq(from=0,to=1,length=alpha)
  }
  if(missing(ncores)){
    if(length(betahat)>50000) ncores=detectCores()
    #Set the number of cores equal to system capacity
  }
  mixcompdist=match.arg(mixcompdist)

  allash=list()
  loglikvector=rep(NA,length(alpha))

  if(ncores==FALSE){
    ##Usual loop without parallel computation
    sink("/dev/null")
    for(i in 1:length(alpha)){
      betahati= betahat/(sebetahat^alpha[i])
      sebetahati= sebetahat^(1-alpha[i])
      beta.ash=ash(betahati, sebetahati, mixcompdist=mixcompdist,df=df,model="EE",...)
      allash[[i]]=beta.ash
      loglikvector[i]=calc_loglik(beta.ash,betahat,sebetahat,df,alpha=alpha[i])
    }
    sink()
  } else{
    ##Performing parallel computation
    cl <- makePSOCKcluster(ncores)#This number corresponding to number of workers
    registerDoParallel(cl)
    allash=foreach(i=1:length(alpha)) %dopar% {
      sink("/dev/null")
      betahati= betahat/(sebetahat^alpha[i])
      sebetahati= sebetahat^(1-alpha[i])
      beta.ash=ashr::ash(betahati, sebetahati, mixcompdist=mixcompdist,df=df,model="EE",...)
      sink()
      beta.ash #computation result stored in allash
    }
    stopCluster(cl)
    for(i in 1:length(alpha)){
      loglikvector[i]=calc_loglik(allash[[i]],betahat,sebetahat,df,alpha=alpha[i])
    }
  }

  modelindex=which.max(loglikvector)
  beta.ash= allash[[modelindex]]
  model=alpha[modelindex]
  if(model==0){
    model="EE"
  } else if(model==1){
    model="ET"
  } else{
    model=model
  }
  beta.ash[["model"]]=model
  return(list(bestash = beta.ash, model=model,loglikevector = loglikvector,allash = allash))
}


#' @title Non-zero-mode Adaptive Shrinkage function
#'
#' @description This is a wrapper function that takes a grid value of
#'     \eqn{mu} and then consider the model \eqn{betahat_j - mu ~
#'     g()},and eqn{beta_j ~ N(0,sebetahat^2) or student t
#'     distribution}. \eqn{mu} should be a grid of mu,this wrapper
#'     function would select the best \eqn{mu} and reports the ash
#'     item based on that \eqn{mu}.
#'
#' @seealso \code{\link{ash}} the main function that this wrapper
#'     function is calling
#' @details All other inputs are exactly the same as the main function
#'     ash, and would pass to the main function to evaluate the
#'     likelihood.
#'
#' @param betahat a p vector of estimates
#' @param sebetahat a p vector of corresponding standard errors
#' @param mixcompdist distribution of components in mixture (
#'     "uniform","halfuniform" or "normal"), the default value would
#'     be "uniform"
#' @param df appropriate degrees of freedom for (t) distribution of
#'     betahat/sebetahat, default is NULL(Gaussian)
#' @param mu Could be a vector of grid values for mu. that this
#'     wrapper would select based on likelihood principle. Could also
#'     be a positive integer greater or equal to 4, then mu number of
#'     grid values would be generated from
#'     [-abs(mean(betahat)),2*abs(mean(betahat)], equally spaced.
#' @param ncores Whether to use parallel computing, defaults to FALSE,
#'     user could specify number of cores they would like to
#'     use. Further, if user does not specify and
#'     length(betahat)>50000, then the function would perform parallel
#'     computation using number of CPU cores on the current host.
#' @param ... Further arguments to be passed to \code{\link{ash}}.
#'
#' @return ashn returns a list of objects \item{beta.ash}{the best
#'     fitted ash object} \item{BestMode}{the best fitted mode, note
#'     that all models are fitted with betahat subtracting the
#'     corresponding mode} \item{loglikvector}{the vector of
#'     loglikelihood of various models} \item{allash}{the fitted ash
#'     of various models}
#'
#' @export
#' @examples
#' beta = c(rep(0,100),rnorm(100))+0.2
#' sebetahat = abs(rnorm(200,0,1))
#' betahat = rnorm(200,beta,sebetahat)
#' beta.ashn = ashn(betahat, sebetahat,mu=20)
#' #beta.ashn4 = ashn(betahat, sebetahat,mu=20,ncores=4)
#' print(beta.ashn[[1]])  #best ash object
#' print(beta.ashn[[2]])  #corresponding mode (0 or some other values)
#' print(beta.ashn[[3]])  #log-likelihood for all models
#'
#'
ashn=function(betahat,sebetahat,
              mixcompdist = c("uniform","halfuniform","normal","+uniform","-uniform"),
			  df=NULL, mu=20,ncores=FALSE,
              ...){
  if(length(mu)==1){
    mutemp=seq(from=(mean(betahat)-3*stats::sd(betahat)/sqrt(length(betahat))),to=(mean(betahat)+3*stats::sd(betahat)/sqrt(length(betahat))),length=mu)
    mu=c(0,mutemp)#include the default 0 mode for comparison
  }

  if(missing(ncores)){
    if(length(betahat)>50000) ncores=detectCores()
    #Set the number of cores equal to system capacity
  }

  mixcompdist = match.arg(mixcompdist)

  allash=list()
  loglikvector=rep(NA,length(mu))

  if(ncores==FALSE){
    ##Usual loop without parallel computation
    sink("/dev/null")
    for(i in 1:length(mu)){
      betahati= betahat-mu[i]
      sebetahati= sebetahat
      beta.ash=ash(betahati, sebetahati, mixcompdist=mixcompdist,df=df,model="EE",...)
      allash[[i]]=beta.ash
      loglikvector[i]=calc_loglik(beta.ash,betahati,sebetahati,df)
    }
    sink()
  } else{
    ##Performing parallel computation
    cl <- makePSOCKcluster(ncores)#This number corresponding to number of workers
    registerDoParallel(cl)
    allash=foreach(i=1:length(mu)) %dopar% {
      sink("/dev/null")
      betahati= betahat-mu[i]
      sebetahati= sebetahat
      beta.ash=ashr::ash(betahati, sebetahati, mixcompdist=mixcompdist,df=df,model="EE",...)
      sink()
      beta.ash #computation result stored in allash
    }
    stopCluster(cl)
    for(i in 1:length(mu)){
      loglikvector[i]=calc_loglik(allash[[i]],betahati,sebetahati,df)
    }
  }

  modelindex=which.max(loglikvector)
  beta.ash=allash[[modelindex]]
  BestMode=mu[modelindex]

  return(list(bestash = beta.ash, BestMode=BestMode,loglikevector = loglikvector,allash = allash))
}


