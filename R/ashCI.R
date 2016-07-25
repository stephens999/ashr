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
ashci = function (a,level=0.95,betaindex,lfsrcriteria=0.05,tol=1e-5, maxcounts=100,shrinkingcoefficient=0.9,trace=FALSE,ncores=FALSE){
  data = a$data
  if(is.null(data)){stop("ash object has to have data returned to compute CIs; use outputlevel 2 or more when running ash")}
 
  options(warn=-1)
  if(missing(betaindex)){
    betaindex = which(get_lfsr(a)<=lfsrcriteria)
    #betaindex[is.na(betaindex)]=FALSE #Some lfsrs would have NA
  }
  
  PosteriorMean = get_pm(a)
  PosteriorSD = get_psd(a)
  ZeroProb = get_lfdr(a)
  NegativeProb = get_np(a)
  PositiveProb = get_pp(a)
  
  m=get_fitted_g(a)
  percentage=1
  
  if( class(m) != "normalmix" && class(m) != "unimix" ){stop(paste("Invalid class",class(m)))}

  CImatrix=matrix(NA,nrow=length(PosteriorMean),ncol=7)
  CImatrix[,3]=PosteriorMean
  colnames(CImatrix)=c("Index(Location)","lfsr","Posterior Mean",(1-level)/2,(1+level)/2,"Fitted CDF(lower) ","Fitted CDF(upper)")
  
  if(missing(trace)){
    if(length(betaindex)>=1000){
      trace=TRUE  #component-wise computation takes more time
    }else {trace=FALSE}
  }
  
  if(missing(ncores)){
    if(length(betaindex)>1000) ncores=detectCores()
  }
  if(ncores==FALSE){
    ## Proceed with sequential computation
    if(trace==TRUE){
      cat("Computation time will be linear w.r.t sample size, progress will be printed to the screen \n")
      tic()
    }
    for(i in betaindex){
      #Starting at Posterior Mean, step out until exceed required level
      #The search will go from the PosteriorMean to the upper or lower point
      #Note: this assumes the Posterior Mean is in the CI... !
      data_i = extract_data(data,i)
      if(is.nan(PosteriorSD[i])){
        CImatrix[i,4]=NA;
      } else if(PosteriorSD[i]==0){ #special case where posterior is (approximately) point mass
        CImatrix[i,4]=PosteriorMean[i];
      } else if(NegativeProb[i]<(1-level)/2 & (ZeroProb[i]+NegativeProb[i])> (1-level)/2){
        CImatrix[i,4]=0;
      } else {
        gap = ifelse(PosteriorSD[i]<1e-5,1e-5,PosteriorSD[i])
        lower = PosteriorMean[i]-gap
        while(cdf_post(m,lower,data_i) > (1-level)/2){
          lower = lower-gap
        }
        CImatrix[i,4]=stats::optimize(f=ci.lower,interval=c(lower,PosteriorMean[i]),m=m,data=data_i,level=level,tol=tol)$minimum
      }
     
      
      if(is.nan(PosteriorSD[i])){
        CImatrix[i,5]=NA;
      } else if(PosteriorSD[i]==0){ #special case where posterior is point mass
        CImatrix[i,5]=PosteriorMean[i];
      } else if(PositiveProb[i] < ((1-level)/2) & (ZeroProb[i]+PositiveProb[i])> (1-level)/2){
        CImatrix[i,5]=0;
      } else {
        gap = ifelse(PosteriorSD[i]<1e-5,1e-5,PosteriorSD[i])
        upper = PosteriorMean[i]+gap
        while(cdf_post(m,upper,data_i) < 1 - (1-level)/2){
          upper = upper+gap
        }
        CImatrix[i,5]=stats::optimize(f=ci.upper,interval=c(PosteriorMean[i],upper),m=m,data=data_i,level=level,tol=tol)$minimum
      }
      
      CImatrix[i,6]=cdf_post(m,CImatrix[i,4],data_i)
      CImatrix[i,7]=cdf_post(m,CImatrix[i,5],data_i)
      
      if(trace==TRUE & percentage <=100){
        currentpercentage=round(i*100/length(betaindex))
        if(currentpercentage == percentage){
          cat("Current computation progress", percentage,"%, seconds ")
          toc()
          percentage = percentage + 1}
      }
    }
    CImatrix[,4] = CImatrix[,4] * data$s_orig^data$alpha
    CImatrix[,5] = CImatrix[,5] * data$s_orig^data$alpha #correct CIs for the fact they are CIs for beta/s^alpha
  } else{
    ## Proceed with parallel computation
    #if(trace==TRUE){
    #cat("Computation time would be linear w.r.t sample size, parallel computation progress would be printed to the screen \n")
    #tic()
    #}
    percentagevector=rep(1,ncores)
    cl <- makePSOCKcluster(ncores)#This number corresponding to number of workers
    registerDoParallel(cl)
    CImatrix[,4:7]=foreach(j=1:length(betaindex), .combine='rbind') %dopar% {
      i = betaindex[j]
      lower = PosteriorMean[i]-PosteriorSD[i]
      data_i = extract_data(data,i)
      while(cdf_post(m,lower,data_i) > (1-level)/2){
        lower = lower-PosteriorSD[i]
      }
      upper = PosteriorMean[i]+PosteriorSD[i]
      while(cdf_post(m,upper,data_i) < 1 - (1-level)/2){
        upper = upper+PosteriorSD[i]
      }
      
      
      #Calculating the lower bound
      CIentryl=stats::optimize(f=ci.lower,interval=c(lower,PosteriorMean[i]),m=m,data=data_i,level=level,tol=tol)$minimum
      cdfl=cdf_post(m, CIentryl,data_i)
      #If the actual cdf deviates by too much, refine the optimization search interval
      #currently we set maximum value of execution to maxcounts to avoid dead loop
      counts=0
      intervallength=PosteriorMean[i]-lower
      while(abs(cdfl-(1-level)/2)>(10*tol) & counts<maxcounts){
        intervallength= intervallength*shrinkingcoefficient
        CIentryl=stats::optimize(f=ci.lower,interval=c(PosteriorMean[i]-intervallength,
                                                       PosteriorMean[i]),m=m,data=data_i,level=level, tol=tol)$minimum
        cdfl=cdf_post(m, CIentryl,data_i)
        counts=counts+1
      }
      
      #Calculating the upper bound
      CIentryu=stats::optimize(f=ci.upper,interval=c(PosteriorMean[i],upper),m=m,data=data_i,level=level,
                               tol=tol)$minimum
      cdfu=cdf_post(m, CIentryu,data_i)
      #If the actual cdf deviates by too much, refine the optimization search interval
      #currently we set maximum value of execution to maxcounts to avoid dead loop
      counts=0
      intervallength=upper-PosteriorMean[i]
      while(abs(cdfu-(1+level)/2)>(10*tol) & counts<maxcounts){
        intervallength= intervallength*shrinkingcoefficient
        CIentryu=stats::optimize(f=ci.upper,interval=c(PosteriorMean[i],PosteriorMean[i]+intervallength),
                                 m=m,data=data_i,level=level,
                                 tol=tol)$minimum
        cdfu=cdf_post(m, CIentryu,data_i)
        counts=counts+1
      }
      
      #sending the result back to master
      c(CIentryl, CIentryu,cdfl,cdfu)
    }
    stopCluster(cl)
  }
  # if(model=="ET"){
  #   CImatrix=CImatrix*sebetahat.orig
  # }
  CImatrix[,1]=1:length(get_lfsr(a))
  CImatrix[,2]=get_lfsr(a)
  CImatrix=signif(CImatrix,digits=round(1-log(tol)/log(10)))
  #CImatrix[,6:7]=round(CImatrix[,6:7],5)
  return(CImatrix)
}

ci.lower=function(z,m,data,level){
  tailprob=cdf_post(m,z,data)
  return((tailprob-(1-level)/2)^2)
}

ci.upper=function(z,m,data,level){
  tailprob=1-cdf_post(m,z,data)
  return((tailprob-(1-level)/2)^2)
}



