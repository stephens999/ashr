

############################### METHODS FOR unimix class ###########################

#' @title Constructor for unimix class
#'
#' @description Creates an object of class unimix (finite mixture of
#'     univariate uniforms)
#'
#' @details None
#'
#' @param pi vector of mixture proportions
#' @param a vector of left hand ends of uniforms
#' @param b vector of right hand ends of uniforms
#'
#' @return an object of class unimix
#'
#' @export
#'
#' @examples unimix(c(0.5,0.5),c(0,0),c(1,2))
unimix = function(pi,a,b){
  structure(data.frame(pi,a,b),class="unimix")
}

#' @export
comp_cdf.unimix = function(m,y,lower.tail=TRUE){
  vapply(y,stats::punif,m$a,min=m$a,max=m$b,lower.tail)
}

comp_sd.unimix = function(m){
  (m$b-m$a)/sqrt(12)
}

#' @export
comp_mean.unimix = function(m){
  (m$a+m$b)/2
}



compdens.unimix = function(m,y,log=FALSE){
  k=ncomp(m)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(stats::dunif(d, m$a, m$b, log),nrow=k))
}

#density of convolution of each component of a unif mixture with s*t_nu() at x
# x an n-vector
#return a k by n matrix
compdens_conv.unimix = function(m,data){
   if(is.null(data$v)){
    compdens= t(stats::pnorm(outer(data$x,m$a,FUN="-")/data$s)-
                  stats::pnorm(outer(data$x,m$b,FUN="-")/data$s))/(m$b-m$a)
    compdens[m$a==m$b,]=t(stats::dnorm(outer(data$x,m$a,FUN="-")/data$s)/data$s)[m$a==m$b,]
  }
  else{
    compdens= t(stats::pt(outer(data$x,m$a,FUN="-")/data$s,df=data$v)-
                  stats::pt(outer(data$x,m$b,FUN="-")/data$s,df=data$v))/(m$b-m$a)
    compdens[m$a==m$b,]=t(stats::dt(outer(data$x,m$a,FUN="-")/data$s,df=data$v)/data$s)[m$a==m$b,]
  }
  return(compdens)
}

#log density of convolution of each component of a unif mixture with s*t_nu() at x
# x an n-vector
#return a k by n matrix
log_compdens_conv.unimix = function(m,data){
  b = pmax(m$b,m$a) #ensure a<b
  a = pmin(m$b,m$a)
  if(is.null(data$v)){
    lpa = stats::pnorm(outer(data$x,a,FUN="-")/data$s,log=TRUE)
    lpb = stats::pnorm(outer(data$x,b,FUN="-")/data$s,log=TRUE)
    
    lcompdens= t( lpa + log(1-exp(lpb-lpa)) ) - log(b-a)
    lcompdens[a==b,]=t(stats::dnorm(outer(data$x,a,FUN="-")/data$s,log=TRUE) - log(data$s))[a==b,]
  }
  else{
    lpa = stats::pt(outer(data$x,a,FUN="-")/data$s,df=data$v,log=TRUE)
    lpb = stats::pt(outer(data$x,b,FUN="-")/data$s,df=data$v,log=TRUE)
    lcompdens= t( lpa + log(1-exp(lpb-lpa)) ) -log(b-a)
    lcompdens[a==b,]=
      t(stats::dt(outer(data$x,a,FUN="-")/data$s,df=data$v,log=TRUE) - log(data$s))[a==b,]
  }
  return(lcompdens)
}


#' @export
compcdf_post.unimix=function(m,c,data){
  k = length(m$pi)
  x = data$x
  s = data$s
  v = data$v
  n=length(x)
  tmp = matrix(1,nrow=k,ncol=n)
  tmp[m$a > c,] = 0
  subset = m$a<=c & m$b>c # subset of components (1..k) with nontrivial cdf
  if(sum(subset)>0){
    if(is.null(v)){
      pna = stats::pnorm(outer(x,m$a[subset],FUN="-")/s)
      pnc = stats::pnorm(outer(x,rep(c,sum(subset)),FUN="-")/s)
      pnb = stats::pnorm(outer(x,m$b[subset],FUN="-")/s)
    }else{
      pna = stats::pt(outer(x,m$a[subset],FUN="-")/s, df=v)
      pnc = stats::pt(outer(x,rep(c,sum(subset)),FUN="-")/s, df=v)
      pnb = stats::pt(outer(x,m$b[subset],FUN="-")/s, df=v)
    }
    tmp[subset,] = t((pnc-pna)/(pnb-pna))
  }
  subset = (m$a == m$b) #subset of components with trivial cdf
  tmp[subset,]= rep(m$a[subset] <= c,n)
  #Occasionally we would encounter issue such that in some entries pna[i,j]=pnb[i,j]=pnc[i,j]=0 or pna=pnb=pnc=1
  #Those are the observations with significant betahat(small sebetahat), resulting in pnorm() return 1 or 0
  #due to the thin tail property of normal distribution.(or t-distribution, although less likely to occur)
  #Then R would be dividing 0 by 0, resulting  in NA values
  #In practice, those observations would have 0 probability of belonging to those "problematic" components
  #Thus any sensible value in [0,1] would not matter much, as they are highly unlikely to come from those
  #components in posterior distribution.
  #Here we simply assign the "naive" value as as (c-a)/(b-a)
  #As the component pdf is rather smaller over the region.
  tmpnaive=matrix(rep((c-m$a)/(m$b-m$a),length(x)),nrow=k,ncol=n)
  tmp[is.nan(tmp)]= tmpnaive[is.nan(tmp)]
  tmp
}

#' @title my_etruncnorm
#' @description Compute expectation of truncated normal.
#'
#' @param a Left limit of distribution.
#' @param b Right limit of distribution.
#' @param mean The mean of the untruncated normal.
#' @param sd The standard deviation of the untruncated normal.
#' @export
my_etruncnorm= function(a,b,mean=0,sd=1){
  alpha = (a-mean)/sd
  beta =  (b-mean)/sd
  #Flip the onese where both are positive, as the computations are more stable
  #when both negative
  flip = (alpha>0 & beta>0)
  flip[is.na(flip)]=FALSE #deal with NAs
  alpha[flip]= -alpha[flip]
  beta[flip]=-beta[flip]
  
  #Fix a bug of quoting the truncnorm package
  #E(X|a<X<b)=a when a==b as a natural result
  #while etruncnorm would simply return NaN,causing PosteriorMean also NaN
  tmp1=etruncnorm(alpha,beta,0,1)
  isequal=(alpha==beta)
  tmp1[isequal]=alpha[isequal]
  
  tmp= (-1)^flip * (mean+sd*tmp1)
  
  max_alphabeta = ifelse(alpha<beta, beta,alpha)
  max_ab = ifelse(alpha<beta,b,a)
  toobig = max_alphabeta<(-30)
  toobig[is.na(toobig)]=FALSE
  tmp[toobig] = max_ab[toobig]
  tmp
}

#' More about the truncated normal
#' @inheritParams my_etruncnorm
#' @export
my_e2truncnorm= function(a,b,mean=0,sd=1){
  alpha = (a-mean)/sd
  beta =  (b-mean)/sd
  #Flip the onese where both are positive, as the computations are more stable
  #when both negative
  flip = (alpha>0 & beta>0)
  flip[is.na(flip)]=FALSE #deal with NAs
  alpha[flip]= -alpha[flip]
  beta[flip]=-beta[flip]
  
  #Fix a bug of quoting the truncnorm package
  #E(X|a<X<b)=a when a==b as a natural result
  #while etruncnorm would simply return NaN,causing PosteriorMean also NaN
  tmp1=etruncnorm(alpha,beta,0,1)
  isequal=(alpha==beta)
  tmp1[isequal]=alpha[isequal]
  tmp= (-1)^flip * (mean+sd*tmp1)
  # for the variance
  # error report in vtruncnorm
  # vtruncnorm(10,-10,0,1)
  # vtruncnorm(-10,10,0,1)
  # vtruncnorm(3,-3,0,1)
  # vtruncnorm(-3,3,0,1)
  # I am not sure smaller one should be put in the first or not
  # vtruncnorm(-7,-8,0,1)
  # vtruncnorm(-8,-7,0,1)
  # vtruncnorm(7,8,0,1)
  # vtruncnorm(8,7,0,1)
  # vtruncnorm(-8,-9,0,1)
  # vtruncnorm(-9,-10,0,1)
  # maybe we should try ourselves according to some result
  # https://people.sc.fsu.edu/~jburkardt/presentations/truncated_normal.pdf
  # tmpvar = vtruncnorm(alpha,beta,0,1)
  tmpvar = my_vtruncnorm(ifelse(alpha<beta,alpha,beta),ifelse(alpha<beta,beta,alpha),0,1)
  # for the second moment
  tmp2 = tmp^2 + tmpvar*sd^2
  isequal=(a==b)
  tmp2[isequal]=(a[isequal])^2
  
  # if the truncate value is too big
  max_alphabeta = ifelse(alpha<beta, beta,alpha)
  max_ab = ifelse(alpha<beta,b,a)
  # I think here, 8 or 7 is enough for this case. try the following:
  toobig = max_alphabeta<(-20)
  toobig[is.na(toobig)]=FALSE
  tmp2[toobig] = (max_ab[toobig])^2
  tmp2
}
#pnorm also have some problems......
#stats::pnorm(7);stats::pnorm(8)
#stats::pnorm(-7);stats::pnorm(-8)
# but it is fine since we flip the sign if a and b are all positive in function my_e2truncnorm
#' This version is better than trunvnorm package
#'
#' @inheritParams my_etruncnorm
#' @export
#'
my_vtruncnorm = function(a,b,mean = 0, sd = 1){
  a = ifelse(a== -Inf, -1e5,a)
  b = ifelse(b== Inf, 1e5, b)
  alpha = (a-mean)/sd
  beta =  (b-mean)/sd
  frac1 = (beta*stats::dnorm(beta,0,1) - alpha*stats::dnorm(alpha,0,1)) / (stats::pnorm(beta,0,1)-stats::pnorm(alpha,0,1) )
  frac2 = (stats::dnorm(beta,0,1) - stats::dnorm(alpha,0,1)) / (stats::pnorm(beta,0,1)-stats::pnorm(alpha,0,1) )
  truncnormvar = sd^2 * (1 - frac1 - frac2^2)
  return(truncnormvar)
}

#note that with uniform prior, posterior is truncated normal, so
#this is computed using formula for mean of truncated normal
#' @export
comp_postmean.unimix = function(m,data){
  #   k= ncomp(m)
  #   n=length(betahat)
  #   a = matrix(m$a,nrow=n,ncol=k,byrow=TRUE)
  #   b = matrix(m$b,nrow=n,ncol=k,byrow=TRUE)
  #   matrix(etruncnorm(a,b,betahat,sebetahat),nrow=k,byrow=TRUE)
  #note: etruncnorm is more stable for a and b negative than positive
  #so maybe use this, and standardize to make the whole more stable.
  x=data$x
  s=data$s
  v=data$v
  alpha = outer(-x, m$a,FUN="+")/s
  beta = outer(-x, m$b, FUN="+")/s
  if(is.null(v)){
    tmp = x + s*my_etruncnorm(alpha,beta,0,1)
  }else{
    tmp = x + s*my_etrunct(alpha,beta,v)
  }
  ismissing = is.na(x) | is.na(s)
  tmp[ismissing,]= (m$a+m$b)/2
  t(tmp)
}

# as for posterior mean, but compute posterior mean squared value
#' @export
comp_postmean2.unimix = function(m,data){
  x=data$x
  s=data$s
  v=data$v
  alpha = outer(-x, m$a,FUN="+")/s
  beta = outer(-x, m$b, FUN="+")/s
  if(is.null(v)){
    tmp = x^2 + 2*x*s*my_etruncnorm(alpha,beta,0,1) + s^2*my_e2truncnorm(alpha,beta,0,1)
  }else{
    tmp = x^2 + 2*x*s*my_etrunct(alpha,beta,v) + s^2*my_e2trunct(alpha,beta,v)
  }
  ismissing = is.na(x) | is.na(s)
  tmp[ismissing,]= (m$b^2+m$a*m$b+m$a^2)/3
  t(tmp)
}

#not yet implemented!
#just returns 0s for now
comp_postsd.unimix = function(m,data){
  k= ncomp(m)
  n=length(data$x)
  return(matrix(NA,nrow=k,ncol=n))
  #  return(sqrt(comp_postmean2(m,betahat,sebetahat,v)-comp_postmean(m,betahat,sebetahat,v)^2))
}

#' @title my_etrunct
#' @description Compute expectation of truncated t, the result is from
#'     the paper
#'     "Moments of truncated Student-t distribution by Hea-Jung Kim"
#'
#' @param a left limit of distribution
#' @param b right limit of distribution
#' @param v degree of freedom of error distribution
#' @export
my_etrunct= function(a,b,v){
  mult=ifelse(a>0 & b>0,-1,1) # calculations more stable for negative
  #a and b, but in this case need to multiply result by -1
  aa = ifelse(a>0 & b>0, -b, a)
  bb = ifelse(a>0 & b>0, -a, b)
  
  A = v+aa^2
  B = v+bb^2
  F_a = stats::pt(aa,df=v)
  F_b = stats::pt(bb,df=v)
  lG=lgamma((v-1)/2)+(v/2)*log(v)-log(2*(F_b-F_a))-lgamma(v/2)-lgamma(1/2)
  G=exp(lG)
  ABpart=(A^(-(v-1)/2)-B^(-(v-1)/2))
  tmp = ifelse(G==Inf & ABpart==0, my_etruncnorm(aa,bb),G*ABpart) #deal with extreme cases using normal
  tmp = ifelse(aa==bb,aa,tmp) #deal with case a=b
  return(mult*tmp)
}
# this is my_e2trunct is wrong function
# my_e2trunct= function(a,b,v){
#   A = v+a^2
#   B = v+b^2
#   F_a = stats::pt(a,df=v)
#   F_b = stats::pt(b,df=v)
#   lG=lgamma((v-1)/2)+(v/2)*log(v)-log(2*(F_b-F_a))-lgamma(v/2)-lgamma(1/2)
#   G=exp(lG)
#   ABpart = (a*A^(-(v-1)/2)-b*B^(-(v-1)/2))
#   # for the second moment
#   EY2 = v/(v-2) + G * ABpart
#   #EY2 = ifelse(G==Inf & ABpart==0, my_e2truncnorm(a,b),EY2) #deal with extreme cases using normal
#   EY2 = ifelse(G==Inf & abs(ABpart)<1e-20, my_e2truncnorm(a,b),EY2)
#   # maybe we also need to deal with extreme case using normal, so I add a truncate normal later
#   return(ifelse(a==b,a^2,EY2)) #deal with extreme case a=b
# }

# we propose a new my_e2trunct function depending on library("hypergeo")
# D_const is a function that used in my_e2trunct, but not been used any more
# D_const = function(A,B,v){
#   f_1 = (A)/(beta(1/2,v/2) * sqrt(v+A^2))
#   f_2 = hypergeo(1/2,1-v/2,3/2,A^2/(v+A^2))
#   part_1 = f_1 * f_2
#   f_1 = (B)/(beta(1/2,v/2) * sqrt(v+B^2))
#   f_2 = hypergeo(1/2,1-v/2,3/2,B^2/(v+B^2))
#   part_2 = f_1 * f_2
#   output = part_1 - part_2
#   return(output)
# }
# this function can be use as any moment calculation, but here we just use it as second moment.
#' @title my_etrunct
#' @description Compute second moment of the truncated t. The result is from the paper "Moments of truncated t and F distributions" by Saralees Nadarajah and Samuel Kotz.
#' @param n is moment, the default is 2 which means second moment
#' @param a left limit of distribution
#' @param b right limit of distribution
#' @param v degree of freedom of error distribution
#' @export
my_e2trunct = function(a,b,v,n=2){
  if(v<=2){warning("my_e2trunct known to be unstable for degrees of freedom v<=2; proceed with caution")}
  aa = ifelse(a>0 & b>0, -b, a)   # calculations more stable for negative
  bb = ifelse(a>0 & b>0, -a, b)  #a and b, because involve cdf(a)-cdf(b)
  a=aa; b=bb; # so switch them where necessary
  
  # deal with infinity case
  a = ifelse(a< (-1e6),-1e6,a)
  b = ifelse(b> (1e6), 1e6, b)
  
  # deal with extreme case, use normal
  # this is just for second moment
  if(v >= 1e5){
    return(my_e2truncnorm(a,b))
  }
  B = a
  A = b
  # D = D_const(A,B,v)
  D = stats::pt(A,df = v) - stats::pt(B,df = v)
  f_1 = 1/((n+1)*sqrt(v)*beta(v/2,1/2)*D)
  f_2 = A^(n+1) * hypergeo::hypergeo((1+v)/2,(1+n)/2,(3+n)/2,-A^2/v) - B^(n+1) * hypergeo::hypergeo((1+v)/2,(1+n)/2,(3+n)/2,-B^2/v)
  output = f_1 * f_2
  
  # deal with same limits case
  output= ifelse(a==b,a^n,output)
  #output = ifelse(Im(output)==0,Re(output),output)
  return(Re(output))
}
