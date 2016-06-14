
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
  print(object$fitted.g)
  print(utils::tail(object$fit$loglik,1),digits=10)
  print(object$fit$converged)
}

#' @title Create data from from ash object
#'
#' @description Creates data frame for easy plotting of results etc
#'
#' @details Returns a data frame with named columnns
#' @param x the fitted ash object
#' @param row.names NULL or a character vector giving the row names for the data frame. Missing values are not allowed.
#' @param optional not used, included for consistency as an S3 generic/method
#' @param ... not used, included for consistency as an S3
#'     generic/method.
#'
#' @examples
#' a = ash(rnorm(100,0,2),1)
#' head(as.data.frame(a))
#' 
#' @export
as.data.frame.ash=function(x,row.names=NULL,optional=FALSE,...){
  if(is.null(x$lfsr)){stop("Can't make data.frame from ash object unless outputlevel>=1.5")}
  if(is.null(row.names)){row.names=1:length(x$lfsr)}
  df = data.frame(row.names=row.names)

  if(!is.null(x$data)){
    df$betahat = x$data$betahat
    df$sebetahat = x$data$sebetahat
  }

  df$lfsr = x$lfsr
  df$svalue = x$svalue
  df$PosteriorMean = x$PosteriorMean
  df$PosteriorSD = x$PosteriorSD

  df$lfdr = x$lfdr
  df$qvalue = x$qvalue

  return(df)
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
  print(x$fitted.g)
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

#compute the predictive density of an observation
#given the fitted ash object a and the vector se of standard errors
#not implemented yet
#todo
predictive=function(a,se){

}


#' @title Get fitted loglikelihood for ash object
#'
#' @description Return the log-likelihood of the data under the fitted
#'     distribution
#'
#' @param a the fitted ash object
#'
#' @details None
#'
#' @export
#'
#'
get_loglik = function(a){
  return(a$loglik)
}

#' @title Get pi0 estimate for ash object
#'
#' @description Return estimate of the null proportion, pi0
#'
#' @param a the fitted ash object
#'
#' @details Extracts the estimate of the null proportion, pi0, from
#'     the object a
#'
#' @export
#'
get_pi0 = function(a){
  null.comp = comp_sd(a$fitted.g)==0
  return(sum(a$fitted.g$pi[null.comp]))
}

#' @title Compute loglikelihood for data from ash fit
#'
#' @description Return the log-likelihood of the data betahat, with
#'     standard errors betahatsd, for a given g() prior on beta, or an
#'     ash object containing that
#'
#' @param g the fitted g, or an ash object containing g
#' @param betahat the data
#' @param betahatsd the observed standard errors
#' @param df appropriate degrees of freedom for (t) distribution of
#'     betahat/sebetahat
#' @param model indicates whether you want the likelihood under the EE
#'     or ET model
#' @param alpha a scalar performing transformation on betahat and
#'     sebetahat, such that the model is \eqn{\beta_j / s_j^alpha ~
#'     g()},and eqn{betahat_j / s_j^alpha ~
#'     N(0,(sebetahat^(1-alpha))^2) or student t distribution}. When
#'     \eqn{alpha=0} we have the EE model, when \eqn{alpha=1}, we have
#'     the ET model. \eqn{alpha} should be in between 0 and 1,
#'     inclusively. Default is alpha=0.
#' @details See example in CompareBetahatvsZscoreAnalysis.rmd
#'
#' @export
#'
calc_loglik = function(g,betahat,betahatsd,df,model=c("EE","ET"),alpha=0){
  if(missing(df)){
    stop("error: must supply df for calc_loglik; supply df=NULL for normal likelihood")
  }

  if(missing(alpha)){
    model = match.arg(model)
    if(class(g)=="ash" && g$model != model){
      warning("Model used to fit ash does not match model used to compute loglik! Probably you have made a mistake!")
    }
    if(model=="ET"){ alpha=1
    } else {alpha=0}
  }

  if(class(g)=="ash"){g = g$fitted.g} #extract g object from ash object if ash object passed


  return(loglik_conv(g,betahat/(betahatsd^alpha),betahatsd^(1-alpha),df)-alpha*sum(log(betahatsd)))
}


#' @title Compute loglikelihood for data under null that all beta are 0
#'
#' @description Return the log-likelihood of the data betahat, with
#'     standard errors betahatsd, under the null that beta==0
#'
#' @param betahat the data
#' @param betahatsd the observed standard errors
#' @param df appropriate degrees of freedom for (t) distribution of
#'     betahat/sebetahat
#' @param model indicates whether you want the likelihood under the EE
#'     or ET model
#' @param alpha a scalar performing transformation on betahat and
#'     sebetahat, such that the model is \eqn{\beta_j / s_j^alpha ~
#'     g()},and eqn{betahat_j / s_j^alpha ~
#'     N(0,(sebetahat^(1-alpha))^2) or student t distribution}. When
#'     \eqn{alpha=0} we have the EE model, when \eqn{alpha=1}, we have
#'     the ET model. \eqn{alpha} should be in between 0 and 1,
#'     inclusively. Default is alpha=0.
#'
#' @export
#'
#'
calc_null_loglik = function(betahat,betahatsd,df,model=c("EE","ET"),alpha=0){
  if(is.null(df)){return(sum(-log(betahatsd) + stats::dnorm(betahat/betahatsd,log=TRUE)))}
  else{return(sum(-log(betahatsd) + stats::dt(betahat/betahatsd,df,log=TRUE)))}
  #g=unimix(1,0,0)
  #return(calc_loglik(g,betahat,betahatsd,df,model,alpha))
}

#' @title Compute loglikelihood ratio for data from ash fit
#'
#' @description Return the log-likelihood ratio of the data betahat,
#'     with standard errors betahatsd, for a given g() prior on beta,
#'     or an ash object containing that, vs the null that g() is point
#'     mass on 0
#'
#' @param g the fitted g, or an ash object containing g
#' @param betahat the data
#' @param betahatsd the observed standard errors
#' @param df appropriate degrees of freedom for (t) distribution of
#'     betahat/sebetahat
#' @param model indicates whether you want the likelihood under the EE
#'     or ET model
#' @param alpha a scalar performing transformation on betahat and
#'     sebetahat, such that the model is \eqn{\beta_j / s_j^alpha ~
#'     g()},and eqn{betahat_j / s_j^alpha ~
#'     N(0,(sebetahat^(1-alpha))^2) or student t distribution}. When
#'     \eqn{alpha=0} we have the EE model, when \eqn{alpha=1}, we have
#'     the ET model. \eqn{alpha} should be in between 0 and 1,
#'     inclusively. Default is alpha=0.
#' @details See example in CompareBetahatvsZscoreAnalysis.rmd
#'
#' @export
calc_logLR = function(g,betahat,betahatsd,df,model=c("EE","ET"),alpha=0){
  return(calc_loglik(g,betahat,betahatsd,df,model,alpha) - calc_null_loglik(betahat,betahatsd,df,model,alpha))
}

#' @title Compute vector of loglikelihood for data from ash fit
#'
#' @description Return the vector of log-likelihoods of the data
#'     betahat, with standard errors betahatsd, for a given g() prior
#'     on beta, or an ash object containing that
#'
#' @param g the fitted g, or an ash object containing g
#' @param betahat the data
#' @param betahatsd the observed standard errors
#' @param df appropriate degrees of freedom for (t) distribution of
#'     betahat/sebetahat
#' @param model indicates whether you want the likelihood under the EE
#'     or ET model
#' @param alpha a scalar performing transformation on betahat and
#'     sebetahat, such that the model is \eqn{\beta_j / s_j^alpha ~
#'     g()},and eqn{betahat_j / s_j^alpha ~
#'     N(0,(sebetahat^(1-alpha))^2) or student t distribution}. When
#'     \eqn{alpha=0} we have the EE model, when \eqn{alpha=1}, we have
#'     the ET model. \eqn{alpha} should be in between 0 and 1,
#'     inclusively. Default is alpha=0.
#' @details See example in CompareBetahatvsZscoreAnalysis.rmd
#'
#' @export
#'
calc_vloglik = function(g,betahat,betahatsd,df,model=c("EE","ET"),alpha=0){
  if(missing(df)){
    stop("error: must supply df for calc_loglik; supply df=NULL for normal likelihood")
  }

  if(missing(alpha)){
    model = match.arg(model)
    if(class(g)=="ash" && g$model != model){
      warning("Model used to fit ash does not match model used to compute loglik! Probably you have made a mistake!")
    }
    if(model=="ET"){ alpha=1
    } else {alpha=0}
  }

  if(class(g)=="ash"){g = g$fitted.g} #extract g object from ash object if ash object passed


  return(log(
    dens_conv(g,betahat/(betahatsd^alpha),betahatsd^(1-alpha),df)-
      alpha*sum(log(betahatsd))))
}


#' @title Compute vector of loglikelihood for data under null that all
#'     beta are 0
#'
#' @description Return the vector of log-likelihoods of the data
#'     betahat, with standard errors betahatsd, under the null that
#'     beta==0
#'
#' @param betahat the data
#' @param betahatsd the observed standard errors
#' @param df appropriate degrees of freedom for (t) distribution of
#'     betahat/sebetahat
#' @param model indicates whether you want the likelihood under the EE
#'     or ET model
#' @param alpha a scalar performing transformation on betahat and
#'     sebetahat, such that the model is \eqn{\beta_j / s_j^alpha ~
#'     g()},and eqn{betahat_j / s_j^alpha ~
#'     N(0,(sebetahat^(1-alpha))^2) or student t distribution}. When
#'     \eqn{alpha=0} we have the EE model, when \eqn{alpha=1}, we have
#'     the ET model. \eqn{alpha} should be in between 0 and 1,
#'     inclusively. Default is alpha=0.
#'
#' @export
#'
#'
calc_null_vloglik = function(betahat,betahatsd,df,model=c("EE","ET"),alpha=0){
  if(is.null(df)){return(-log(betahatsd) + stats::dnorm(betahat/betahatsd,log=TRUE))}
  else{return(-log(betahatsd) + stats::dt(betahat/betahatsd,df,log=TRUE))}
  #g=unimix(1,0,0)
  #return(calc_loglik(g,betahat,betahatsd,df,model,alpha))
}

#' @title Compute vector of loglikelihood ratio for data from ash fit
#'
#' @description Return the vector of log-likelihood ratios of the data
#'     betahat, with standard errors betahatsd, for a given g() prior
#'     on beta, or an ash object containing that, vs the null that g()
#'     is point mass on 0
#'
#' @param g the fitted g, or an ash object containing g
#' @param betahat the data
#' @param betahatsd the observed standard errors
#' @param df appropriate degrees of freedom for (t) distribution of
#'     betahat/sebetahat
#' @param model indicates whether you want the likelihood under the EE
#'     or ET model
#' @param alpha a scalar performing transformation on betahat and
#'     sebetahat, such that the model is \eqn{\beta_j / s_j^alpha ~
#'     g()},and eqn{betahat_j / s_j^alpha ~
#'     N(0,(sebetahat^(1-alpha))^2) or student t distribution}. When
#'     \eqn{alpha=0} we have the EE model, when \eqn{alpha=1}, we have
#'     the ET model. \eqn{alpha} should be in between 0 and 1,
#'     inclusively. Default is alpha=0.
#' @details See example in CompareBetahatvsZscoreAnalysis.rmd
#'
#' @export
calc_vlogLR = function(g,betahat,betahatsd,df,model=c("EE","ET"),alpha=0){
  return(calc_vloglik(g,betahat,betahatsd,df,model,alpha) - calc_null_vloglik(betahat,betahatsd,df,model,alpha))
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
  list(x=x,y=dens(a$fitted.g,x))
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
  return(list(x=x,y=mixcdf(a$fitted.g,x,lower.tail)))
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



#' @title Credible Interval Computation for the ash object
#'
#' @description Given the ash object return by the main function ash,
#'     this function computes the corresponding credible interval of
#'     the mixture model.
#'
#' @details Uses default optimization function and perform
#'     component-wise credible interval computation. The computation
#'     cost is linear of the length of betahat.
#'
#' @param a the fitted ash object
#' @param betahat the values of beta used for ash
#' @param sebetahat the values of betahat used
#' @param df the degrees of freedom
#' @param model the model used.
#' @param level the level for the credible interval, (default=0.95)
#' @param betaindex a vector consisting of locations of betahat where
#'     you would like to compute the credible interval
#' @param lfsrcriteria a scalar, in which the function would
#'     autoselect betahats based on lfsr value smaller than
#'     lfsrcriteria when index is not supplied. Setting it to 1 would
#'     compute credible interval for all observations.
#' @param tol the desired accuracy, default value is 1e-5.
#' @param maxcounts a positive integer used as the maximum iterations
#'     to carry additional optimization on a refined interval,when you
#'     see Actual cdf deviates from (1-level)/2,(1+level)/2 by too
#'     much, try to increase this number. (default is 100)
#' @param shrinkingcoefficient A positive real number smaller than 1
#'     (default is 0.9), used to shrink to search interval for lower
#'     and upper bound of the credible interval
#' @param trace a logical variable denoting whether some of the
#'     intermediate results of iterations should be displayed to the
#'     user. Default is FALSE.
#' @param ncores Whether to use parallel computing, defaults to FALSE,
#'     user could specify number of cores they would like to
#'     use. Further, if user does not specify and we have over 1000
#'     entries to compute, then the function would perform parallel
#'     computation using number of CPU cores on the current host.
#' @return A matrix, with first column being the location,second
#'     column being lfsr, 3rd column being posterior mean,4th and 5th
#'     column being the lower bound and upper bound for the credible
#'     interval.
#'
#'
#' @export
#' @examples
#' beta = c(rep(0,20),rnorm(20))
#' sebetahat = abs(rnorm(40,0,1))
#' betahat = rnorm(40,beta,sebetahat)
#' beta.ash = ash(betahat, sebetahat)
#'
#' CImatrix=ashci(beta.ash,betahat,sebetahat,level=0.95)
#' print(CImatrix)
#' print(CImatrix[order(CImatrix[,2]),]) # Sorted according to the lfsr
#'
#' CImatrix1=ashci(beta.ash,betahat,sebetahat,level=0.95,betaindex=c(1,2,5))
#' CImatrix2=ashci(beta.ash,betahat,sebetahat,level=0.95,lfsrcriteria=0.1)
#' #CImatrix3=ashci(beta.ash,betahat,sebetahat,level=0.95, betaindex=c(1:length(beta)),ncores=4)
#' print(CImatrix1)
#' print(CImatrix2)
#' #print(CImatrix3)
#'
#' ##A larger example
#' #beta = c(rep(0,1000),rnorm(1000))
#' #sebetahat = abs(rnorm(2000,0,1))
#' #betahat = rnorm(2000,beta,sebetahat)
#' #beta.ash = ash(betahat, sebetahat)
#' #CImatrix4 = ashci(beta.ash,betahat,sebetahat,level=0.95, betaindex=c(1:length(beta)),ncores=4)
#todo/issue=> all set!
#1.Q:Could do parallel computing to reduce the computation time
#1.A:Done by doParallel
#
#2.Q:Optimization function does not work well for discountinous function, even
#if it's a monotone one. Current remedy is to set a more conservative value for
#the searching interval from the mixture
#2.A:Done by shrinking searching interval using while loop
#
ashci = function (a, betahat=NULL, sebetahat=NULL,df=NULL,model=c("EE","ET"),level=0.95,betaindex,lfsrcriteria=0.05,tol=1e-5, maxcounts=100,shrinkingcoefficient=0.9,trace=FALSE,ncores=FALSE){
  if(missing(betahat)){
      betahat= a$data$betahat
  }
  if(missing(sebetahat)){
      sebetahat = a$data$sebetahat
  }
  if(missing(df)){
      df = a$data$df
  }

  options(warn=-1)
  if(missing(betaindex)){
    betaindex =(a$lfsr<=lfsrcriteria)
    betaindex[is.na(betaindex)]=FALSE #Some lfsrs would have NA
  }

  PosteriorMean=a$PosteriorMean[betaindex]
  PosteriorSD = a$PosteriorSD[betaindex]
  ZeroProb = a$ZeroProb[betaindex]
  NegativeProb = a$NegativeProb[betaindex]
  PositiveProb = a$PositiveProb[betaindex]
  x=betahat[betaindex]
  s=sebetahat[betaindex]
  m=a$fitted.g
  model=match.arg(model)
  percentage=1

  if( class(m) != "normalmix" && class(m) != "unimix" ){stop(paste("Invalid class",class(m)))}
  if(model=="ET"){ #for ET model, standardize
    x=x/s
    PosteriorMean=PosteriorMean/s
    sebetahat.orig=s
    s=rep(1,length(x))
  }
  if(is.null(df)){
    errorspan=stats::qnorm(level)
  } else{
    errorspan=stats::qt(level,df)
  }

  CImatrix=matrix(NA,nrow=length(x),ncol=7)
  CImatrix[,3]=PosteriorMean
  colnames(CImatrix)=c("Index(Location)","lfsr","Posterior Mean",(1-level)/2,(1+level)/2,"Fitted CDF(lower) ","Fitted CDF(upper)")

  if(missing(trace)){
    if(length(x)>=1000){
      trace=TRUE  #component-wise computation takes more time
    }else {trace=FALSE}
  }

  if(missing(ncores)){
    if(length(x)>1000) ncores=detectCores()
  }
  if(ncores==FALSE){
    ## Proceed with sequential computation
    if(trace==TRUE){
      cat("Computation time would be linear w.r.t sample size, progress would be printed to the screen \n")
      tic()
    }
    for(i in 1:length(x)){
    #Starting at Posterior Mean, step out until exceed required level
      #The search will go from the PosteriorMean to the upper or lower point
      #Note: this assumes the Posterior Mean is in the CI... !
      lower = PosteriorMean[i]-PosteriorSD[i]
      while(cdf_post(m,lower,x[i],s[i],df) > (1-level)/2){
          lower = lower-PosteriorSD[i]
      }
      upper = PosteriorMean[i]+PosteriorSD[i]
      while(cdf_post(m,upper,x[i],s[i],df) < 1 - (1-level)/2){
        upper = upper+PosteriorSD[i]
      }


      #Calculating the lower bound
      #First check if lower bound is 0
      if(NegativeProb[i]<(1-level)/2 & (ZeroProb[i]+NegativeProb[i])> (1-level)/2){
          CImatrix[i,4]=0;
          CImatrix[i,6]=cdf_post(m,CImatrix[i,4],x[i],s[i],df)
      } else {
        CImatrix[i,4]=stats::optimize(f=ci.lower,interval=c(lower,PosteriorMean[i]),m=m,x=x[i],s=s[i],level=level,
                             df=df, tol=tol)$minimum
        CImatrix[i,6]=cdf_post(m,CImatrix[i,4],x[i],s[i],df)
      }
      #Calculating the upper bound
      #First check if upper bound is 0
      if(PositiveProb[i] < ((1-level)/2) & (ZeroProb[i]+PositiveProb[i])> (1-level)/2){
        CImatrix[i,5]=0;
        CImatrix[i,7]=cdf_post(m,CImatrix[i,5],x[i],s[i],df)
      } else {
        CImatrix[i,5]=stats::optimize(f=ci.upper,interval=c(PosteriorMean[i],upper),m=m,x=x[i],s=s[i],level=level,
                             df=df, tol=tol)$minimum
        CImatrix[i,7]=cdf_post(m,CImatrix[i,5],x[i],s[i],df)
      }
      if(trace==TRUE & percentage <=100){
        currentpercentage=round(i*100/length(x))
        if(currentpercentage == percentage){
          cat("Current computation progress", percentage,"%, seconds ")
          toc()
          percentage = percentage + 1}
      }
    }
  } else{
    ## Proceed with parallel computation
    #if(trace==TRUE){
    #cat("Computation time would be linear w.r.t sample size, parallel computation progress would be printed to the screen \n")
    #tic()
    #}
    percentagevector=rep(1,ncores)
    cl <- makePSOCKcluster(ncores)#This number corresponding to number of workers
    registerDoParallel(cl)
    CImatrix[,4:7]=foreach(i=1:length(x), .combine='rbind') %dopar% {
      lower = PosteriorMean[i]-PosteriorSD[i]
      while(cdf_post(m,lower,x[i],s[i],df) > (1-level)/2){
        lower = lower-PosteriorSD[i]
      }
      upper = PosteriorMean[i]+PosteriorSD[i]
      while(cdf_post(m,upper,x[i],s[i],df) < 1 - (1-level)/2){
        upper = upper+PosteriorSD[i]
      }


      #Calculating the lower bound
      CIentryl=stats::optimize(f=ci.lower,interval=c(lower,PosteriorMean[i]),m=m,x=x[i],s=s[i],level=level,
                        df=df, tol=tol)$minimum
      cdfl=cdf_post(m, CIentryl,x[i],s[i],df)
      #If the actual cdf deviates by too much, refine the optimization search interval
      #currently we set maximum value of execution to maxcounts to avoid dead loop
      counts=0
      intervallength=PosteriorMean[i]-lower
      while(abs(cdfl-(1-level)/2)>(10*tol) & counts<maxcounts){
        intervallength= intervallength*shrinkingcoefficient
        CIentryl=stats::optimize(f=ci.lower,interval=c(PosteriorMean[i]-intervallength,
                                                       PosteriorMean[i]),m=m,x=x[i],s=s[i],level=level,df=df, tol=tol)$minimum
        cdfl=cdf_post(m, CIentryl,x[i],s[i],df)
        counts=counts+1
      }

      #Calculating the upper bound
      CIentryu=stats::optimize(f=ci.upper,interval=c(PosteriorMean[i],upper),m=m,x=x[i],s=s[i],level=level,
                        df=df, tol=tol)$minimum
      cdfu=cdf_post(m, CIentryu,x[i],s[i],df)
      #If the actual cdf deviates by too much, refine the optimization search interval
      #currently we set maximum value of execution to maxcounts to avoid dead loop
      counts=0
      intervallength=upper-PosteriorMean[i]
      while(abs(cdfu-(1+level)/2)>(10*tol) & counts<maxcounts){
        intervallength= intervallength*shrinkingcoefficient
        CIentryu=stats::optimize(f=ci.upper,interval=c(PosteriorMean[i],PosteriorMean[i]+intervallength),
                          m=m,x=x[i],s=s[i],level=level,
                          df=df, tol=tol)$minimum
        cdfu=cdf_post(m, CIentryu,x[i],s[i],df)
        counts=counts+1
      }

      #sending the result back to master
      c(CIentryl, CIentryu,cdfl,cdfu)
    }
    stopCluster(cl)
  }
  if(model=="ET"){
    CImatrix=CImatrix*sebetahat.orig
  }
  numericindex=c(1:length(a$data$betahat))[betaindex]
  CImatrix[,1]=numericindex
  CImatrix[,2]=a$lfsr[betaindex]
  CImatrix=signif(CImatrix,digits=round(1-log(tol)/log(10)))
  #CImatrix[,6:7]=round(CImatrix[,6:7],5)
  return(CImatrix)
}

ci.lower=function(z,m,x,s,level,df){
  tailprob=cdf_post(m,z,x,s,df)
  return((tailprob-(1-level)/2)^2)
}

ci.upper=function(z,m,x,s,level,df){
  tailprob=1-cdf_post(m,z,x,s,df)
  return((tailprob-(1-level)/2)^2)
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


#' @title Multi-Mode Adaptive Shrinkage function
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
#' @return ashm returns a list of objects \item{beta.ash}{the best
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
