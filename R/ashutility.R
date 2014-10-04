
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
#' @param zscores indicates whether ash object was originally fit to z scores 
#' @details None
#' 
#' @export
#' 
#'
loglik.ash = function(a,betahat,betahatsd,zscores=FALSE){
  g=a$fitted.g
  FUN="+"
  if(zscores==TRUE){
    g$sd = sqrt(g$sd^2+1) 
    FUN="*"
  }
  return(loglik_conv(g,betahat, betahatsd,FUN))
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
#' @param levels, the level for the credible interval, (default=0.95)
#' 
#' @return A matrix, with first column being the posterior mean, second and third column being the lower bound and upper bound for the credible interval. 
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
#' 
#' 
#' 
ashci = function (a,level=0.95){
  options(warn=-1)
  x=a$data$betahat
  s=a$data$sebetahat
  m=a$fitted.g
  df=a$df
  if(class(m)=="normalmix"){
    lower=min(x)-qt(level,df=1)*(max(m$sd)+max(abs(s)))
    upper=max(x)+qt(level,df=1)*(max(m$sd)+max(abs(s)))
  } else{
    lower=min(c(m$a,m$b))
    upper=max(c(m$a,m$b))
  }

  CImatrix=matrix(NA,nrow=length(x),ncol=3)	
  colnames(CImatrix)=c("Posterior Mean",(1-level)/2,(1+level)/2)
  CImatrix[,1]=a$PosteriorMean
	
  if( class(a$fitted.g) == "normalmix" | class(a$fitted.g) == "unimix" ){
    for(i in 1:length(x)){
	  CImatrix[i,2]=optim(par=a$PosteriorMean[i],f=ci.lower,m=m,x=x[i],s=s[i],level=level,df=df,method="Brent",lower=lower,upper=upper)$par
	  CImatrix[i,3]=optim(par=a$PosteriorMean[i],f=ci.upper,m=m,x=x[i],s=s[i],level=level,df=df,method="Brent",lower=lower,upper=upper)$par
	}
  } else{stop(paste("Invalid class",class(m)))}
  return(CImatrix)
}

ci.lower=function(z,m,x,s,level,df){
	tailprob=cdf_post(m,z,x,s,df)
	return(abs(tailprob-(1-level)/2))
}

ci.upper=function(z,m,x,s,level,df){
	tailprob=1-cdf_post(m,z,x,s,df)
	return(abs(tailprob-(1-level)/2))
}


