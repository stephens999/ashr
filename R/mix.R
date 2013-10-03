
################################## GENERIC FUNCTIONS ############################
# find matrix of densities at y, for each component of the mixture
# INPUT y is an n-vector
# OUTPUT k by n matrix of densities
compdens = function(x,y,log=FALSE){
  UseMethod("compdens")
}
compdens.default = function(x,y,log=FALSE){
  stop("No such class")
}

#standard deviations
comp_sd = function(m){
  UseMethod("comp_sd")
}
comp_sd.default = function(m){
  stop("method comp_sd not written for this class")
}

#return the overall mean of the mixture
mixmean = function(m){
  UseMethod("mixmean")
}
mixmean.default = function(m){
  sum(m$pi * comp_mean(m))
}

#standard deviations
comp_mean = function(m){
  UseMethod("comp_mean")
}
comp_mean.default = function(m){
  stop("method comp_mean not written for this class")
}

#number of components
ncomp = function(m){
  UseMethod("ncomp")
}
ncomp.default = function(m){
  return(length(m$pi))
}

#return mixture proportions, a generic function
mixprop = function(m){
  UseMethod("mixprop")
}
mixprop.default = function(m){
  m$pi
}

#find cdf at y, a generic function
mixcdf = function(x,y,lower.tail=TRUE){
  UseMethod("mixcdf")
}
mixcdf.default = function(x,y,lower.tail=TRUE){
  x$pi %*% comp_cdf(x,y,lower.tail)
}

#find cdf for each component, a generic function
comp_cdf = function(x,y,lower.tail=TRUE){
  UseMethod("comp_cdf")
}
comp_cdf.default = function(x,y,lower.tail=TRUE){
  stop("comp_cdf not implemented for this class")
}


#find density at y, a generic function
dens = function(x,y){
  UseMethod("dens")
}
dens.default = function(x,y){
  return (x$pi %*% compdens(x, y))
}

#find log likelihood of data in y
LogLik = function(x,y){
  UseMethod("LogLik")
}
LogLik.default = function(x,y){
  sum(dens(x,y))
}

#compute the density of the components of the mixture m
#when convoluted with a normal with standard deviation s
#the density is evaluated at x
#x and s are n-vectors
#m is a mixture with k components
#output is a k by n matrix of densities
compdens_conv = function(m, x, s){
  UseMethod("compdens_conv")
}
compdens_conv.default = function(m,x, s){
  stop("No such class")
}

#compute density of mixture m convoluted with normal of sd (s)
#at locations x
#m is a mixture
#x is an n vector
#s is an n vector or integer
dens_conv = function(m,x,s){
  UseMethod("dens_conv")
}
dens_conv.default = function(m,x,s){
  colSums(m$pi * compdens_conv(m,x,s))
}

#compute the posterior prob that each observation
#came from each component of the mixture m
#output a k by n vector of probabilities
#computed by weighting the component densities by pi
#and then normalizing
comppostprob=function(m,x,s){
 UseMethod("comppostprob") 
}
comppostprob.default = function(m,x,s){
  tmp= (t(m$pi * compdens_conv(m,x,s))/dens_conv(m,x,s))
  ismissing = (is.na(x) | is.na(s))
  tmp[ismissing,]=m$pi
  t(tmp)
}
# evaluate cdf of posterior distribution of beta at c
# m is the prior on beta, a mixture
# c is location of evaluation
# assumption is betahat | beta \sim N(beta,sebetahat)
# m is a mixture with k components
# c a scalar
# betahat, sebetahat are n vectors 
# output is a k by n matrix
compcdf_post=function(m,c,betahat,sebetahat){
  UseMethod("compcdf_post")
}
compcdf_post.default=function(m,c,betahat,sebetahat){
  stop("method compcdf_post not written for this class")
}


cdf_post = function(m,c,betahat,sebetahat){
  UseMethod("cdf_post")
}
cdf_post.default=function(m,c,betahat,sebetahat){
  colSums(comppostprob(m,betahat,sebetahat)*compcdf_post(m,c,betahat,sebetahat))
}

#output posterior mean for beta for prior mixture m,
#given observations betahat, sebetahat
postmean = function(m, betahat,sebetahat){
  UseMethod("postmean")
}
postmean.default = function(m,betahat,sebetahat){
  colSums(comppostprob(m,betahat,sebetahat) * comp_postmean(m,betahat,sebetahat))
}



#output posterior mean for beta for each component of prior mixture m,
#given observations betahat, sebetahat
comp_postmean = function(m, betahat,sebetahat){
  UseMethod("comp_postmean")
}
comp_postmean.default = function(m,betahat,sebetahat){
  stop("method comp_postmean not written for this class")
}


#find nice limits of mixture m for plotting
min_lim = function(m){
  UseMethod("min_lim")
}
min_lim.default=function(m){
  -5
}

max_lim = function(m){
  UseMethod("max_lim")
}
max_lim.default=function(m){
  5
}


#plot density of mixture
plot_dens = function(m,npoints=100,...){
  UseMethod("plot_dens")
}
plot_dens.default = function(m,npoints=100,...){
  x = seq(min_lim(m),max_lim(m),length=npoints)
  plot(x,dens(m,x),type="l",xlab="density",ylab="x",...)
}

plot_post_cdf = function(m,betahat,sebetahat,npoints=100,...){
  UseMethod("plot_post_cdf")
}
plot_post_cdf.default = function(m,betahat,sebetahat,npoints=100,...){
  x = seq(min_lim(m),max_lim(m),length=npoints)
  x_cdf = vapply(x,cdf_post,FUN.VALUE=betahat,m=m,betahat=betahat,sebetahat=sebetahat)
  plot(x,x_cdf,type="l",xlab="x",ylab="cdf",...)
 # for(i in 2:nrow(x_cdf)){
 #   lines(x,x_cdf[i,],col=i)
 # }
}

############################### METHODS FOR normalmix class ###########################

# constructor
normalmix = function(pi,mean,sd){
  structure(data.frame(pi,mean,sd),class="normalmix")
}

comp_sd.normalmix = function(m){
  m$sd
}

comp_mean.normalmix = function(m){
  m$mean
}

compdens.normalmix = function(x,y,log=FALSE){
  k=ncomp(x)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(dnorm(d, x$mean, x$sd, log),nrow=k))  
}

#density of convolution of each component of a normal mixture with N(0,s^2) at x
# x an n-vector at which density is to be evaluated
#return a k by n matrix
#Note that convolution of two normals is normal, so it works that way
compdens_conv.normalmix = function(m, x, s){
  if(length(s)==1){s=rep(s,length(x))}
  sdmat = sqrt(outer(s^2,m$sd^2,FUN="+")) #n by k matrix of standard deviations of convolutions
  return(t(dnorm(outer(x,m$mean,FUN="-")/sdmat)/sdmat))
}


comp_cdf.normalmix = function(x,y,lower.tail=TRUE){
  vapply(y,pnorm,x$mean,x$mean,x$sd,lower.tail)
}

#c is a scalar
#m a mixture with k components
#betahat a vector of n observations
#sebetahat an n vector of standard errors
#return a k by n matrix of the posterior cdf
compcdf_post.normalmix=function(m,c,betahat,sebetahat){
  k = length(m$pi)
  n=length(betahat)
  #compute posterior standard deviation (s1) and posterior mean (m1)
  s1 = sqrt(outer(sebetahat^2,m$sd^2,FUN="*")/outer(sebetahat^2,m$sd^2,FUN="+"))
  ismissing = (is.na(betahat) | is.na(sebetahat))
  s1[ismissing,]=m$sd
  
  m1 = t(comp_postmean(m,betahat,sebetahat))
  t(pnorm(c,mean=m1,sd=s1))
}

#return posterior mean for each component of prior m, given observations betahat and sebetahat
#input, m is a mixture with k components
#betahat, sebetahat are n vectors
#output is a k by n matrix
comp_postmean.normalmix = function(m,betahat,sebetahat){
  tmp=(outer(sebetahat^2,m$mean, FUN="*") + outer(betahat,m$sd^2, FUN="*"))/outer(sebetahat^2,m$sd^2,FUN="+")
  ismissing = (is.na(betahat) | is.na(sebetahat))
  tmp[ismissing,]=m$mean #return prior mean when missing data
  t(tmp)
}



############################### METHODS FOR unimix class ###########################

#constructor
#here the uniform is parameterized in terms of min=mean-sd and min=mean+sd
#(so obviously sd is a misnomer!)
unimix = function(pi,a,b){
  structure(data.frame(pi,a,b),class="unimix")
}

comp_cdf.unimix = function(m,y,lower.tail=TRUE){
  vapply(y,punif,m$a,min=m$a,max=m$b,lower.tail)
}

comp_sd.unimix = function(m){
  (m$b-m$a)/sqrt(12)
}

comp_mean.unimix = function(m){
  (m$a+m$b)/2
}

compdens.unimix = function(x,y,log=FALSE){
  k=ncomp(x)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(dunif(d, x$a, x$b, log),nrow=k))  
}

#density of convolution of each component of a unif mixture with N(0,s) at x
# x an n-vector
#return a k by n matrix
compdens_conv.unimix = function(m, x, s){
  return(t(pnorm(outer(x,m$a,FUN="-")/s)
          -pnorm(outer(x,m$b,FUN="-")/s))/(m$b-m$a))
}


#c is a scalar
#m a mixture with k components
#betahat a vector of n observations
#sebetahat an n vector of standard errors
#return a k by n matrix of the posterior cdf
compcdf_post.unimix=function(m,c,betahat,sebetahat){
  k = length(m$pi)
  n=length(betahat)
  tmp = matrix(1,nrow=k,ncol=n)
  tmp[m$a >= c,] = 0
  subset = m$a<c & m$b>c # subset of components (1..k) with nontrivial cdf
  if(sum(subset)>0){
    pna = pnorm(outer(betahat,m$a[subset],FUN="-")/sebetahat)
    pnc = pnorm(outer(betahat,rep(c,sum(subset)),FUN="-")/sebetahat)
    pnb = pnorm(outer(betahat,m$b[subset],FUN="-")/sebetahat)
    tmp[subset,] = t((pnc-pna)/(pnb-pna))
  }
  tmp
}

my_etruncnorm= function(a,b,mean=0,sd=1){
  alpha = (a-mean)/sd
  beta =  (b-mean)/sd
 #Flip the onese where both are positive, as the computations are more stable
  #when both negative
  flip = (alpha>0 & beta>0)
  flip[is.na(flip)]=FALSE #deal with NAs
  alpha[flip]= -alpha[flip]
  beta[flip]=-beta[flip]
  
  tmp= (-1)^flip * (mean+sd*etruncnorm(alpha,beta,0,1))
  
  max_alphabeta = ifelse(alpha<beta, beta,alpha)
  max_ab = ifelse(alpha<beta,b,a)
  toobig = max_alphabeta<(-30)
  toobig[is.na(toobig)]=FALSE 
  tmp[toobig] = max_ab[toobig]
  tmp
}
  
  
#return posterior mean for each component of prior m, given observations betahat and sebetahat
#input, m is a mixture with k components
#betahat, sebetahat are n vectors
#output is a k by n matrix
#note that with uniform prior, posterior is truncated normal, so
#this is computed using formula for mean of truncated normal 
comp_postmean.unimix = function(m,betahat,sebetahat){
#   k= ncomp(m)
#   n=length(betahat)
#   a = matrix(m$a,nrow=n,ncol=k,byrow=TRUE)
#   b = matrix(m$b,nrow=n,ncol=k,byrow=TRUE)
#   matrix(etruncnorm(a,b,betahat,sebetahat),nrow=k,byrow=TRUE)
  #note: etruncnorm is more stable for a and b negative than positive
  #so maybe use this, and standardize to make the whole more stable.
  
  alpha = outer(-betahat, m$a,FUN="+")/sebetahat
  beta = outer(-betahat, m$b, FUN="+")/sebetahat
  tmp = betahat + sebetahat*my_etruncnorm(alpha,beta,0,1)
  ismissing = is.na(betahat) | is.na(sebetahat)
  tmp[ismissing,]= (m$a+m$b)/2
  t(tmp)
#   t(
#     betahat + sebetahat* 
#       exp(dnorm(alpha,log=TRUE)- pnorm(alpha,log=TRUE))
#    * 
#       (-expm1(dnorm(beta,log=TRUE)-dnorm(alpha,log=TRUE)))
#     /
#       (expm1(pnorm(beta,log=TRUE)-pnorm(alpha,log=TRUE)))
#   )
}


