
################################## ASH UTILITY FUNCTIONS ############################

#' @title Summary method for ash object
#'
#' @description Print summary of fitted ash object
#'
#' @details \code{\link{summary}} prints the fitted mixture, the fitted log likelihood with 10 digits and a flag to indicate convergence
#'
#' @export
#' 
summary.ash=function(a){
  print(a$fitted.g)
  print(tail(a$fit$loglik,1),digits=10)
  print(a$fit$converged)
}

#' @title Print method for ash object
#'
#' @description Print the fitted distribution of beta values in the EB hierarchical model
#'
#' @details None
#' 
#' @export
#' 
print.ash =function(a){
  print(a$fitted.g)
}

#' @title Plot method for ash object
#'
#' @description Plot the density of the underlying fitted distribution
#'
#' @details None
#' 
#' @export
#' 
plot.ash = function(a,xmin,xmax,...){
  x = seq(xmin,xmax,length=1000)
  y = density(a,x)
  plot(y,type="l",...)
}

#compute the predictive density of an observation
#given the fitted ash object a and the vector se of standard errors
#not implemented yet
predictive=function(a,se){
  
}


#' @title Get fitted loglikelihood for ash object
#'
#' @description Return the log-likelihood of the data under the fitted distribution
#'
#' @param a the fitted ash object
#'
#' @details None
#' 
#' @export
#' 
#'
get_loglik = function(a){
  return(tail(a$fit$loglik,1))
}

#' @title Get pi0 estimate for ash object
#'
#' @description Return estimate of the null proportion, pi0
#'
#' @param a the fitted ash object
#'
#' @details Extracts the estimate of the null proportion, pi0, from the object a
#' 
#' @export
#' 
get_pi0 = function(a){
  null.comp = comp_sd(a$fitted.g)==0
  return(sum(a$fitted.g$pi[null.comp]))
}

#' @title Compute loglikelihood for data from ash fit
#'
#' @description Return the log-likelihood of the data betahat, with standard errors betahatsd, 
#' under the fitted distribution in the ash object. 
#' 
#'
#' @param a the fitted ash object
#' @param betahat the data
#' @param betahatsd the observed standard errors
#' @param model: indicates whether you want the likelihood under the EE or ES model 
#' @details See example in CompareBetahatvsZscoreAnalysis.rmd
#' 
#' @export
#' 
#'
calc_loglik = function(a,betahat,betahatsd,df,model=c("EE","ES")){
  if(missing(df)){
    stop("error: must supply df for calc_loglik")
  }
  model = match.arg(model) 
  if(a$model != model){
    warning("Model used to fit ash does not match model used to compute loglik! Probably you have made a mistake!")
  }
  g=a$fitted.g
  if(model=="ES"){
      return(loglik_conv(g,betahat/betahatsd,1,df)-sum(log(betahatsd)))
  } else {
      return(loglik_conv(g,betahat,betahatsd,df))
  }
}

#' @title Compute loglikelihood for data usign the prior g (usually from an ash fit)
#'
#' @description Return the log-likelihood of the data betahat, with standard errors betahatsd, 
#' under the fitted distribution in the ash object. 
#' 
#'
#' @param g the prior for effects or standardized effects
#' @param betahat the data
#' @param betahatsd the observed standard errors
#' @param model: indicates whether you want the likelihood under the EE or ES model 
#' @details See example in CompareBetahatvsZscoreAnalysis.rmd
#' 
#' @export
#' 
#'
calc_gloglik = function(g,betahat,betahatsd,df,model=c("EE","ES")){
  if(missing(df)){
    stop("error: must supply df for calc_loglik")
  }
  model = match.arg(model) 
  if(model=="ES"){
    return(loglik_conv(g,betahat/betahatsd,1,df)-sum(log(betahatsd)))
  } else {
    return(loglik_conv(g,betahat,betahatsd,df))
  }
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
density.ash=function(a,x){list(x=x,y=dens(a$fitted.g,x))}

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
#' 
#'
cdf.ash=function(a,x,lower.tail=TRUE){
  return(list(x=x,y=mixcdf(a$fitted.g,x,lower.tail)))
}



#' @title Credible Interval Computation for the ash object
#'
#' @description Given the ash object return by the main function ash, this function computes the corresponding credible interval of the mixture model.
#'
#' @details Uses default optimization function and perform component-wise credible interval computation. The computation cost is linear of the length of betahat.
#'
#' @param a the fitted ash object 
#' @param betaindex a vector consisting of locations of betahat where you would like to compute the credible interval
#' @param lfsrcriteria a scalar, in which the function would autoselect betahats based on lfsr value smaller than lfsrcriteria when index is not supplied. Setting it to 1 would compute credible interval for all observations.
#' @param levels the level for the credible interval, (default=0.95)
#' @param trace a logical variable denoting whether some of the intermediate results of iterations should be displayed to the user. Default is FALSE.
#' @param tol the desired accuracy, default value is 1e-5.
#' 
#' @return A matrix, with first column being the location,second column being lfsr, 3rd column being posterior mean,4th and 5th column being the lower bound and upper bound for the credible interval. 
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
#' CImatrix2=ashci(beta.ash,level=0.95,lfsrcriteria=0.1)
#' 
#' 
#todo/issue
#1.Could do parallel computing to reduce the computation time
#
#2.Optimization function does not work well for discountinous function, even
#if it's a monotone one. Current remedy is to set a more conservative value for 
#the searching interval from the mixture
#
ashci = function (a,level=0.95,betaindex,lfsrcriteria=0.05,tol=1e-5,trace=FALSE){
  #options(warn=-1)
  if(missing(betaindex)){
  	betaindex =(a$lfsr<=lfsrcriteria)
  	betaindex[is.na(betaindex)]=FALSE #Some lfsrs would have NA
  }  
  x=a$data$betahat[betaindex]
  s=a$data$sebetahat[betaindex]
  PosteriorMean=a$PosteriorMean[betaindex]
  m=a$fitted.g
  df=a$df
  model=a$model
  percentage=1
  
  if(model=="ES"){ #for ES model, standardize
  	x=x/s
	PosteriorMean=PosteriorMean/s
  	sebetahat.orig=s
  	s=rep(1,length(x))
  }

  if(missing(trace)){
    if(length(x)>=1000){
      trace=TRUE  #component-wise computation takes more time
    }else {trace=FALSE}
  }
  if(trace==TRUE){
  	cat("Computation time would be linear w.r.t sample size, progress would be printed to the screen \n")
  	tic()
  }
  
  if(is.null(df)){
  	errorspan=qnorm(level)
  } else{
  	errorspan=qt(level,df)
  }
  
  CImatrix=matrix(NA,nrow=length(x),ncol=5)
  CImatrix[,3]=PosteriorMean
  colnames(CImatrix)=c("Index(Location)","lfsr","Posterior Mean",(1-level)/2,(1+level)/2)
  
  if( class(m) == "normalmix" | class(m) == "unimix" ){
    for(i in 1:length(x)){
      #Now the search interval is better restricted, avoiding the crash of optimize() due to discontinuity of cdf_post
      #The discontinuity is due to the pointmass component of the mixture
      #cumpi=cumsum(comppostprob(m,x[i],s[i],df))-level
      #maxposition=min(which(cumpi>0))
      if(class(m)=="normalmix"){
      	#maxsd=m$sd[maxposition]
      	#lower=PosteriorMean[i]-errorspan*maxsd
      	#upper=PosteriorMean[i]+errorspan*maxsd
      	lower=-Inf
      	upper=Inf
	  }else{
        #lower=min(c(m$a[1: maxposition],m$b[1: maxposition]))
        #upper=max(c(m$a[1: maxposition],m$b[1: maxposition])) 	
        lower=min(c(m$a,m$b))
        upper=max(c(m$a,m$b))
	  }
      
      CImatrix[i,4]=optimize(f=ci.lower,interval=c(lower,PosteriorMean[i]),m=m,x=x[i],s=s[i],level=level,
      df=df, tol=tol)$minimum
	  
	  CImatrix[i,5]=optimize(f=ci.upper,interval=c(PosteriorMean[i],upper),m=m,x=x[i],s=s[i],level=level,
	  df=df, tol=tol)$minimum
	  
	  #CImatrix[i,2]=optim(par=a$PosteriorMean[i],f=ci.lower,m=m,x=x[i],s=s[i],level=level,
	  #df=df,method="Brent",lower=lower,upper=upper)$par
	  #CImatrix[i,3]=optim(par=a$PosteriorMean[i],f=ci.upper,m=m,x=x[i],s=s[i],level=level,
	  #df=df,method="Brent",lower=lower,upper=upper)$par	  
	  if(trace==TRUE & percentage <=100){
	  	currentpercentage=round(i*100/length(x))
	  	if(currentpercentage == percentage){
	  		cat("Current computation progress", percentage,"%, seconds ")
	  		toc()
	  		percentage = percentage + 1}
	  }	  
	}
  } else{stop(paste("Invalid class",class(m)))}
  
  if(model=="ES"){
    CImatrix=CImatrix*sebetahat.orig
  }
  numericindex=c(1:length(a$data$betahat))[betaindex]
  CImatrix[,1]=numericindex
  CImatrix[,2]=a$lfsr[betaindex]
  #CImatrix=signif(CImatrix,digits=round(1-log(tol)/log(10)))
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





