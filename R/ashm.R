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


