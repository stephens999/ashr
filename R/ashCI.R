#' @title Credible Interval Computation for the ash object
#'
#' @description Given the ash object returned by the main function ash,
#'     this function computes a posterior credible interval (CI) for each observation. The ash object
#'     must include a data component to use this function (which it does by default).
#'
#' @details Uses uniroot to find credible interval, one at a time for each observation. 
#' The computation cost is linear in number of observations.
#'
#' @param a the fitted ash object
#' @param level the level for the credible interval, (default=0.95)
#' @param betaindex a vector consisting of locations of betahat where
#'     you would like to compute the credible interval
#' @param lfsr_threshold a scalar, if specified then computes CIs only for observations
#' more significant than that threshold. 
#' @param tol passed to uniroot; indicates desired accuracy.
#' @param trace a logical variable denoting whether some of the
#'     intermediate results of iterations should be displayed to the
#'     user. Default is FALSE.

#' @return A matrix, with 2 columns, ith row giving CI for ith observation
#'
#'
#' @export
#' @examples
#' beta = c(rep(0,20),rnorm(20))
#' sebetahat = abs(rnorm(40,0,1))
#' betahat = rnorm(40,beta,sebetahat)
#' beta.ash = ash(betahat, sebetahat)
#'
#' CImatrix=ashci(beta.ash,level=0.95)
#'
#' CImatrix1=ashci(beta.ash,level=0.95,betaindex=c(1,2,5))
#' CImatrix2=ashci(beta.ash,level=0.95,lfsr_threshold=0.1)
ashci = function (a,level=0.95,betaindex,lfsr_threshold=1,tol=1e-3,trace=FALSE){
  data = a$data
  if(is.null(data)){stop("ash object has to have data returned to compute CIs; use outputlevel 2 or more when running ash")}
 
  # options(warn=-1)
  if(missing(betaindex)){
    betaindex = which(get_lfsr(a)<=lfsr_threshold)
    #betaindex[is.na(betaindex)]=FALSE #Some lfsrs would have NA
  }
  
  PosteriorMean = get_pm(a)
  PosteriorSD = get_psd(a)
  ZeroProb = get_lfdr(a)
  NegativeProb = get_np(a)
  PositiveProb = get_pp(a)
  
  m=get_fitted_g(a)
  percentage=1
  
  if( class(m) != "normalmix" && class(m) != "unimix" && class(m) != "tnormalmix"){stop(paste("Invalid class",class(m)))}

  CImatrix=matrix(NA,nrow=length(PosteriorMean),ncol=2)
  colnames(CImatrix)=c((1-level)/2,(1+level)/2)
  #c("Fitted CDF(lower) ","Fitted CDF(upper)")
  
  if(missing(trace)){
    if(length(betaindex)>=1000){
      trace=TRUE  #component-wise computation takes more time
    }else {trace=FALSE}
  }

  if(trace==TRUE){
      cat("Computation time will be linear w.r.t sample size, progress will be printed to the screen \n")
      tic()
  }
  for(i in betaindex){
    data_i = extract_data(data,i)
    
    if(is.nan(PosteriorSD[i])){
      CImatrix[i,]=c(NA,NA)
    } else if(PosteriorSD[i]==0){ #special case where posterior is (approximately) point mass
      CImatrix[i,]=PosteriorMean[i]
    } else {
      #find lower point (first checking if 0 is it)
      if(NegativeProb[i]<(1-level)/2 & (ZeroProb[i]+NegativeProb[i])> (1-level)/2){
        CImatrix[i,1]=0;
      } else {
        CImatrix[i,1]=stats::uniroot(f=taildiff,interval=c(PosteriorMean[i]-2*PosteriorSD[i],PosteriorMean[i]),extendInt="upX",m=m,data=data_i,target=(1-level)/2,tol=tol)$root
      }
      
      #find upper point (first checking if 0 is it)
      if(PositiveProb[i] < ((1-level)/2) & (ZeroProb[i]+PositiveProb[i])> (1-level)/2){
        CImatrix[i,2]=0;
      } else {
        CImatrix[i,2]=stats::uniroot(f=taildiff,interval=c(PosteriorMean[i],PosteriorMean[i]+2*PosteriorSD[i]),extendInt="upX",m=m,data=data_i,target=(1+level)/2,tol=tol)$root
      }
    }
    if(trace==TRUE & percentage <=100){
      currentpercentage=round(i*100/length(betaindex))
      if(currentpercentage == percentage){
        cat("Current computation progress", percentage,"%, seconds ")
        toc()
        percentage = percentage + 1}
    }
  }
  CImatrix = CImatrix * data$s_orig^data$alpha #correct CIs for the fact they are CIs for beta/s^alpha
  
  CImatrix=signif(CImatrix,digits=round(1-log(tol)/log(10)))
  return(CImatrix)
}


#difference of tailprob from target
taildiff=function(z,m,data,target){
  cdf_post(m,z,data)-target
}




