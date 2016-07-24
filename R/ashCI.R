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
  
  PosteriorMean=a$res$PosteriorMean[betaindex]
  PosteriorSD = a$res$PosteriorSD[betaindex]
  ZeroProb = a$res$lfdr[betaindex]
  NegativeProb = a$res$NegativeProb[betaindex]
  PositiveProb = a$res$PositiveProb[betaindex]
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
  CImatrix[,2]=a$res$lfsr[betaindex]
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



