
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

#' @export
ashn2 = function(betahat,...){
  test.op = function(c){return(-ash(betahat=betahat,mode=c,output="loglik",...)$loglik)}
  opt = optimize(test.op,interval=c(min(betahat),max(betahat)))
  ash(betahat=betahat,mode=opt$minimum,...)
}

ash.nzm = function(...){
  cc = sys.call()
  test.op = function(c){eval(call_modify(cc,newargs = list(mode=c)))$loglik}
  test.op(2)
   # opt = optimize(test.op,interval=c(min(betahat),max(betahat)))
#  ash(betahat=betahat,mode=opt$minimum,...)
}
