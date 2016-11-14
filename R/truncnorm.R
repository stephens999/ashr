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
  #E(X|a<X<b)=a when a==b is a natural result
  #while etruncnorm would simply return NaN,causing PosteriorMean also NaN
  # ZMZ:  when a and b are both negative and far from 0, etruncnorm can't compute
  # the mean and variance. Also we should deal with 0/0 situation caused by sd = 0.
  tmp1=etruncnorm(alpha,beta,0,1)
  isequal=is.equal(alpha,beta)
  
  tmp1[isequal]=alpha[isequal]
  
  tmp= (-1)^flip * (mean+sd*tmp1)
  
  max_alphabeta = ifelse(alpha<beta, beta,alpha)
  max_ab = ifelse(alpha<beta,b,a)
  toobig = max_alphabeta<(-30)
  toobig[is.na(toobig)]=FALSE
  tmp[toobig] = max_ab[toobig]

  NAentry = is.na(tmp)
  sdd = sd
  sdd[NAentry] = 0
  tmp = modify.sd0(tmp,a,b,mean,sdd)
  
  tmp
}

#tests for equality, with NA defined to be FALSE
is.equal = function(a,b){
  isequal = (a==b)
  isequal[is.na(isequal)]=FALSE
  return(isequal)
}

#elementwisely deal with entries with sd equals to 0
modify.sd0 = function(temp,a,b,mean,sd){
  dimen = dim(temp)
  if(is.null(dimen)){
    dimen = c(length(temp),1)
  }
  A = matrix(a,dimen[1],dimen[2])
  B = matrix(b,dimen[1],dimen[2])
  M = matrix(mean,dimen[1],dimen[2])
  SD = matrix(sd,dimen[1],dimen[2])
  sdzero = which(sd == 0)
  temp[sdzero] = ifelse(A[sdzero]<=M[sdzero] & B[sdzero]>=M[sdzero],M[sdzero],ifelse(A[sdzero] > M[sdzero],A[sdzero],B[sdzero]))
  temp
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
  isequal=is.equal(alpha,beta)
 
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
  isequal=is.equal(a,b)
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
  
  # turn all nan and negative into 0
  nan.index = is.na(truncnormvar) | truncnormvar < 0
  truncnormvar[nan.index] = 0
  return(truncnormvar)
}
