
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
#' @param levels the level for the credible interval, (default=0.95)
#' @param betaindex a vector consisting of locations of betahat where you would like to compute the credible interval
#' @param lfsrcriteria a scalar, in which the function would autoselect betahats based on lfsr value smaller than lfsrcriteria when index is not supplied. Setting it to 1 would compute credible interval for all observations.
#' @param tol the desired accuracy, default value is 1e-5.
#' @param maxcounts a positive integer used as the maximum iterations to carry additional optimization on a refined interval,when you see Actual cdf deviates from (1-level)/2,(1+level)/2 by too much, try to increase this number. (default is 100) 
#' @param shrinkingcoefficient A positive real number smaller than 1 (default is 0.9), used to shrink to search interval for lower and upper bound of the credible interval
#' @param trace a logical variable denoting whether some of the intermediate results of iterations should be displayed to the user. Default is FALSE.
#' @param ncores Whether to use parallel computing, defaults to FALSE, user could specify number of cores they would like to use. Further, if user does not specify and we have over 1000 entries to compute,  then the function would perform parallel computation  using number of CPU cores on the current host.
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
#' print(CImatrix)
#'
#' CImatrix1=ashci(beta.ash,level=0.95,betaindex=c(1,2,5))
#' CImatrix2=ashci(beta.ash,level=0.95,lfsrcriteria=0.1)
#' CImatrix3=ashci(beta.ash,level=0.95, betaindex=c(1:length(beta)),ncores=4)
#' print(CImatrix1)
#' print(CImatrix2)
#' print(CImatrix3)
#' 
#' #A larger example
#' beta = c(rep(0,1000),rnorm(1000))
#' sebetahat = abs(rnorm(2000,0,1))
#' betahat = rnorm(2000,beta,sebetahat)
#' beta.ash = ash(betahat, sebetahat)
#' CImatrix4 = ashci(beta.ash,level=0.95, betaindex=c(1:length(beta)),ncores=4)
#todo/issue=> all set!
#1.Could do parallel computing to reduce the computation time
#1. Done by doParallel
#
#2.Optimization function does not work well for discountinous function, even
#if it's a monotone one. Current remedy is to set a more conservative value for 
#the searching interval from the mixture
#2.Done by shrinking searching interval using while loop
#
ashci = function (a,level=0.95,betaindex,lfsrcriteria=0.05,tol=1e-5, maxcounts=100,shrinkingcoefficient=0.9,trace=FALSE,ncores=FALSE){
  options(warn=-1)  
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
  
  if( class(m) != "normalmix" && class(m) != "unimix" ){stop(paste("Invalid class",class(m)))}
  if(model=="ES"){ #for ES model, standardize
  	x=x/s
	PosteriorMean=PosteriorMean/s
  	sebetahat.orig=s
  	s=rep(1,length(x))
  }
  if(is.null(df)){
  	errorspan=qnorm(level)
  } else{
  	errorspan=qt(level,df)
  }
  
  CImatrix=matrix(NA,nrow=length(x),ncol=7)
  CImatrix[,3]=PosteriorMean
  colnames(CImatrix)=c("Index(Location)","lfsr","Posterior Mean",(1-level)/2,(1+level)/2,"Actual cdf(lower) ","Actual cdf(upper)")  
  
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
      #Now the search interval is better restricted, avoiding the crash of optimize() due to discontinuity of cdf_post
      #The discontinuity is due to the pointmass component of the mixture
      cumpi=cumsum(comppostprob(m,x[i],s[i],df))-level
      maxposition=min(which(cumpi>0))
      if(class(m)=="normalmix"){
      	maxsd=m$sd[maxposition]
      	lower=PosteriorMean[i]-errorspan*maxsd
      	upper=PosteriorMean[i]+errorspan*maxsd
	  }else{
        lower=min(c(m$a[1: maxposition],m$b[1: maxposition]))
        upper=max(c(m$a[1: maxposition],m$b[1: maxposition])) 	
	  }
	  
	  #Calculating the lower bound
      CImatrix[i,4]=optimize(f=ci.lower,interval=c(lower,PosteriorMean[i]),m=m,x=x[i],s=s[i],level=level,
      df=df, tol=tol)$minimum
	  CImatrix[i,6]=cdf_post(m,CImatrix[i,4],x[i],s[i],df)
	  #If the actual cdf deviates by too much, refine the optimization search interval
	  #currently we set maximum value of execution to maxcounts to avoid dead loop
	  counts=0
	  intervallength=PosteriorMean[i]-lower
	  while(abs(CImatrix[i,6]-(1-level)/2)>(10*tol) & counts<maxcounts){
	  	intervallength= intervallength* shrinkingcoefficient
        CImatrix[i,4]=optimize(f=ci.lower,interval=c(PosteriorMean[i]-intervallength,
        PosteriorMean[i]),m=m,x=x[i],s=s[i],level=level,df=df, tol=tol)$minimum
	    CImatrix[i,6]=cdf_post(m,CImatrix[i,4],x[i],s[i],df)	  
	    counts=counts+1	
	  }
	  
	  #Calculating the upper bound
	  CImatrix[i,5]=optimize(f=ci.upper,interval=c(PosteriorMean[i],upper),m=m,x=x[i],s=s[i],level=level,
	  df=df, tol=tol)$minimum
	  CImatrix[i,7]=cdf_post(m,CImatrix[i,5],x[i],s[i],df)
	  #If the actual cdf deviates by too much, refine the optimization search interval
	  #currently we set maximum value of execution to maxcounts to avoid dead loop
	  counts=0
	  intervallength=upper-PosteriorMean[i]
	  while(abs(CImatrix[i,7]-(1+level)/2)>(10*tol) & counts<maxcounts){
	  	intervallength= intervallength*shrinkingcoefficient
	    CImatrix[i,5]=optimize(f=ci.upper,interval=c(PosteriorMean[i],PosteriorMean[i]+intervallength),m=m,x=x[i],s=s[i],level=level,
	    df=df, tol=tol)$minimum
	    CImatrix[i,7]=cdf_post(m,CImatrix[i,5],x[i],s[i],df)  	
	    counts=counts+1		    
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
      cumpi=cumsum(ashr:::comppostprob(m,x[i],s[i],df))-level
      maxposition=min(which(cumpi>0))
      if(class(m)=="normalmix"){
      	maxsd=m$sd[maxposition]
      	lower=PosteriorMean[i]-errorspan*maxsd
      	upper=PosteriorMean[i]+errorspan*maxsd
	  }else{
        lower=min(c(m$a[1: maxposition],m$b[1: maxposition]))
        upper=max(c(m$a[1: maxposition],m$b[1: maxposition])) 	
	  }
	  
	  #Calculating the lower bound	  
      CIentryl=optimize(f=ashr:::ci.lower,interval=c(lower,PosteriorMean[i]),m=m,x=x[i],s=s[i],level=level,
      df=df, tol=tol)$minimum	  
	  cdfl=cdf_post(m, CIentryl,x[i],s[i],df)
	  #If the actual cdf deviates by too much, refine the optimization search interval
	  #currently we set maximum value of execution to maxcounts to avoid dead loop
	  counts=0
	  intervallength=PosteriorMean[i]-lower	  
	  while(abs(cdfl-(1-level)/2)>(10*tol) & counts<maxcounts){
	  	intervallength= intervallength*shrinkingcoefficient
        CIentryl=optimize(f=ashr:::ci.lower,interval=c(PosteriorMean[i]-intervallength,
        PosteriorMean[i]),m=m,x=x[i],s=s[i],level=level,df=df, tol=tol)$minimum	  
	    cdfl=cdf_post(m, CIentryl,x[i],s[i],df)
	    counts=counts+1	
	  }
	  
	  #Calculating the upper bound
	  CIentryu=optimize(f=ashr:::ci.upper,interval=c(PosteriorMean[i],upper),m=m,x=x[i],s=s[i],level=level,
	  df=df, tol=tol)$minimum
	  cdfu=cdf_post(m, CIentryu,x[i],s[i],df)
	  #If the actual cdf deviates by too much, refine the optimization search interval
	  #currently we set maximum value of execution to maxcounts to avoid dead loop
	  counts=0
	  intervallength=upper-PosteriorMean[i]
	  while(abs(cdfu-(1+level)/2)>(10*tol) & counts<maxcounts){
	  	intervallength= intervallength*shrinkingcoefficient
	    CIentryu=optimize(f=ashr:::ci.upper,interval=c(PosteriorMean[i],PosteriorMean[i]+intervallength),m=m,x=x[i],s=s[i],level=level,
	    df=df, tol=tol)$minimum
	    cdfu=cdf_post(m, CIentryu,x[i],s[i],df)
	    counts=counts+1		    
	  }
	  	  
	  #sending the result back to master
	  c(CIentryl, CIentryu,cdfl,cdfu)
	}
   stopCluster(cl)
  }
  if(model=="ES"){
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
#' @description This is a wrapper function that takes a grid value of \eqn{\alpha} and then consider the model \eqn{\beta_j / s_j^{\alpha} ~ g()},when \eqn{\alpha=0} we have the EE model, when \eqn{\alpha=1}, we have the ES model. \eqn{\alpha} should be in between 0 and 1, inclusively. This wrapper function would select the best \eqn{\alpha} and reports the ash item based on that \eqn{\alpha}.

#'
#' @seealso \code{\link{ash}} the main function that this wrapper function is calling
#' @details All other inputs are exactly the same as the main function ash, and would pass to the main function to evaluate the likelihood.
#'
#' @param betahat  a p vector of estimates 
#' @param sebetahat a p vector of corresponding standard errors
#' @param method specifies how ash is to be run. Can be "shrinkage" (if main aim is shrinkage) or "fdr" (if main aim is to assess fdr or fsr)
#' This is simply a convenient way to specify certain combinations of parameters: "shrinkage" sets pointmass=FALSE and prior="uniform";
#' "fdr" sets pointmass=TRUE and prior="nullbiased".
#' @param mixcompdist distribution of components in mixture ( "uniform","halfuniform" or "normal"), the default value would be "uniform"
#' @param lambda1  multiplicative "inflation factor" for standard errors (like Genomic Control)
#' @param lambda2  additive "inflation factor" for standard errors (like Genomic Control)
##' @param df appropriate degrees of freedom for (t) distribution of betahat/sebetahat, default is NULL(Gaussian)
#' @param nullweight scalar, the weight put on the prior under "nullbiased" specification, see \code{prior}
#' @param nonzeromode logical, indicating whether to use a non-zero unimodal mixture(default is "FALSE")
#' @param alpha Could be a vector of grid values in interval [0,1], that this wrapper would select based on likelihood principle. Could also be a positive integer greater or equal to 2, then alpha number of grid values would be generated from [0,1], equally spaced. The default value is 2 that we compare the EE and ES model.
#' @param ncores Whether to use parallel computing, defaults to FALSE, user could specify number of cores they would like to use. Further, if user does not specify and length(betahat)>50000, then the function would perform parallel computation using number of CPU cores on the current host.
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
#' beta.ashm4 = ashm(betahat, sebetahat,alpha=6,ncores=4)
#' print(beta.ashm[[1]])
#' print(beta.ashm[[2]])
#' 


ashm=function(betahat,sebetahat,method = c("shrink","fdr"), 
               mixcompdist = c("uniform","halfuniform","normal"),
               lambda1=1,lambda2=0,df=NULL,
               nullweight=10,nonzeromode=FALSE, 
               alpha=2,ncores=FALSE){
  if(length(alpha)==1){
  	alpha=seq(from=0,to=1,length=alpha)
  }
  if(missing(ncores)){
  	if(length(betahat)>50000) ncores=detectCores()
  	#Set the number of cores equal to system capacity
  }
  
  if(ncores==FALSE){
  	##Usual loop without parallel computation
    sink("/dev/null")
    allash=list()
    loglikvector=rep(NA,length(alpha))
    for(i in 1:length(alpha)){
      betahati= betahat/(sebetahat^alpha[i])
      sebetahati= sebetahat^(1-alpha[i])	
      beta.ash=ash(betahati, sebetahati, method=method, mixcompdist=mixcompdist, lambda1=lambda1,
      lambda2=lambda2, df=df,nullweight= nullweight,
      nonzeromode= nonzeromode,model="EE")
	  allash[[i]]=beta.ash
	  loglikvector[i]= beta.ash$fit$loglik
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
      beta.ash=ashr::ash(betahati, sebetahati, method=method, mixcompdist=mixcompdist, lambda1=lambda1,
      lambda2=lambda2, df=df,nullweight= nullweight,
      nonzeromode= nonzeromode,model="EE")
	  sink()
	  beta.ash #computation result stored in allash
    }
    stopCluster(cl)
    for(i in 1:length(alpha)){
	  loglikvector[i]= allash[[i]]$fit$loglik
    } 	
  }
  modelindex=which(loglikvector==max(loglikvector))
  beta.ash= allash[[modelindex]]
  model=alpha[modelindex]
  if(model==0){
  	model="EE"
  } else if(model==1){
  	  	model="ES"
  }
  beta.ash[["model"]]=model
  return(list(beta.ash, loglikvector, allash))
}

