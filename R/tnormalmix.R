tnormalmix = function(pi,mean,sd,a,b){
  structure(data.frame(pi,mean,sd,a,b),class="tnormalmix")
}

comp_sd.tnormalmix = function(m){
  #tnormalvar = numeric(length(m$pi))
  #for(k in 1:length(m$pi)){
  #  tnormalvar[k] = my_vtruncnorm(m$a[k],m$b[k],m$mean[k],m$sd[k])
  #}
  tnormalvar = my_vtruncnorm(m$a,m$b,m$mean,m$sd)
  sqrt(tnormalvar)
}

comp_mean.tnormalmix = function(m){
  #tnormalmean = rep(0,length(m$pi))
  #for(k in 1:length(m$pi)){
  #  tnormalmean[k] = my_etruncnorm(m$a[k],m$b[k],m$mean[k],m$sd[k])
  #}
  tnormalmean = my_etruncnorm(m$a,m$b,m$mean,m$sd)
  tnormalmean  
}

compdens.tnormalmix = function(m,y,log=FALSE){
  k = ncomp(m)
  n = length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  ### use package truncnorm
  ##return(matrix(truncnorm::dtruncnorm(d,m$a,m$b,m$mean,m$sd,log),nrow=k))
  ## no cases of b=a
  return(matrix(stats::dnorm(d,m$mean,m$sd))/(stats::pnorm(m$b)-stats::pnorm(m$a)))
}

comp_dens_conv.tnormalmix = function(m,data){
  if(!is_normal(data$lik)){
    stop("Error: truncated normal mixture for non-normal likelihood is not yet implemented")
  }
  if(length(data$s)==1){data$s=rep(data$s,length(data$x))}
  A = sqrt(outer(1/m$sd^2,1/data$s^2,FUN="+"))
  B = 1/sqrt(outer(m$sd^2,data$s^2,FUN="+"))
  C = outer(m$sd,data$s,"/")
  D = stats::pnorm(m$b/m$sd)-stats::pnorm(m$a/m$sd)
  varmat = outer(m$sd^2,data$s^2,FUN="+")
  left = stats::pnorm(A*m$b-t(t(B*C)*data$x))
  right = stats::pnorm(A*m$a-t(t(B*C)*data$x))
  denx = stats::dnorm(matrix(data$x,length(m$sd),length(data$x),byrow=TRUE)/sqrt(varmat))/sqrt(varmat)
  result = ((left-right)*denx)/D
  DD = stats::dnorm(m$b/m$sd)
  lleft = stats::dnorm(A*m$b-t(t(B*C)*data$x))
  result[m$a==m$b,] = (((lleft)*denx/varmat)/DD)[m$a==m$b,]
  result[m$sd==0,] = denx[m$sd==0,]
  return(result)
}

log_comp_dens_conv.tnormalmix = function(m,data) {
  if(!is_normal(data$lik)){
    stop("Error: truncated normal mixture for non-normal likelihood is not yet implemented")
  }
  #### use previous function directly
  #return(log(compdens_conv.tnormalmix(m,x,s,v)))
  if(length(data$s)==1){data$s=rep(data$s,length(data$x))}
  A = sqrt(outer(1/m$sd^2,1/data$s^2,FUN="+"))
  B = 1/sqrt(outer(m$sd^2,data$s^2,FUN="+"))
  C = outer(m$sd,data$s,"/")
  D = stats::pnorm(m$b/m$sd)-stats::pnorm(m$a/m$sd)
  varmat = outer(m$sd^2,data$s^2,FUN="+")
  left = stats::pnorm(A*m$b-t(t(B*C)*data$x))
  right = stats::pnorm(A*m$a-t(t(B*C)*data$x))
  denx = stats::dnorm(t(matrix(data$x,length(data$x),length(m$sd))),0,sqrt(varmat),log=TRUE)
  result = log(left-right)+denx-log(matrix(D,length(m$sd),length(data$x)))
  DD = stats::dnorm(m$b/m$sd)
  lleft = stats::dnorm(A*m$b-t(t(B*C)*data$x))
  result[m$a==m$b,] = (log(lleft/DD)+denx-log(varmat))[m$a==m$b,]
  result[m$sd==0,] = denx[m$sd==0,]
  return(result)
}

#vapply(1:10, pnorm, c(1,2,3), c(1,2,3),c(1,2,3))
comp_cdf.tnormalmix = function(m,y,lower.tail=TRUE){
  k = length(m$pi)
  n=length(y)
  tmp = matrix(1,nrow=k,ncol=n)
  #tmp[m$a > y,] = 0
  subset=outer(m$a,y,'>')
  tmp[subset] = 0
  #subset = m$a<=y & m$b>y
  subset1 = outer(m$a,y,'<=')
  subset2 = outer(m$b,y,'>=')
  subset = subset1 & subset2
  if(sum(subset)>0){
    #pnc=vapply(y,stats::pnorm,m$mean[subset],m$mean[subset],m$sd[subset],lower.tail)
    #pna=vapply(rep(m$a,n),stats::pnorm,m$mean[subset],m$mean[subset],m$sd[subset],lower.tail)
    #pnb=vapply(rep(m$b,n),stats::pnorm,m$mean[subset],m$mean[subset],m$sd[subset],lower.tail)
    YY = matrix(y,k,n,byrow=TRUE)
    MM = matrix(m$mean,k,n)
    SD = matrix(m$sd,k,n)
    pnc = matrix(stats::pnorm(YY[subset],MM[subset],SD[subset]),k,n)
    A=matrix(m$a,k,ncol=n)
    pna=matrix(stats::pnorm(A[subset],MM[subset],SD[subset],lower.tail),k,n)
    B=matrix(m$b,k,ncol=n)
    pnb=matrix(stats::pnorm(B[subset],MM[subset],SD[subset],lower.tail),k,n)
  }
  tmp[subset] = (pnc-pna)/(pnb-pna)
  #subset=(m$a==m$b)
  #tmp[subset,]=rep(m$a[subset]<=c,n)
  tmp
}
  
comp_cdf_post.tnormalmix = function(m,c,data){
  if(!is_normal(data$lik)){
    stop("Error: truncated normal mixture for non-normal likelihood is not yet implemented")
  }  
  k = length(m$pi)
  n=length(data$x)
  tmp = matrix(1,nrow=k,ncol=n)
  tmp[m$a > c,] = 0
  subset = m$a<=c & m$b>=c
  if(sum(subset)>0){
    X=1/(outer(data$s^2,m$sd[subset]^2,FUN="/")+1)
    Y=1/outer(1/data$s^2,1/m$sd[subset]^2,FUN="+")
    A=matrix(m$a[subset],nrow=sum(subset),ncol=n)
    pna=stats::pnorm(t(A),X*data$x+t(t(1-X)*m$mean[subset]),sqrt(Y))
    C=matrix(c,nrow=sum(subset),ncol=n)
    pnc=stats::pnorm(t(C),X*data$x+t(t(1-X)*m$mean[subset]),sqrt(Y))
    B=matrix(m$b[subset],nrow=sum(subset),ncol=n)
    pnb=stats::pnorm(t(B),X*data$x+t(t(1-X)*m$mean[subset]),sqrt(Y))
  }
  tmp[subset,]=t((pnc-pna)/(pnb-pna))
  subset=(m$a==m$b)
  tmp[subset,]=rep(m$a[subset]<=c,n)
  subset=B==C
  tmp[subset] = 1
  ### ZMZ: in case of pnc = pnb, we make it 1 and other
  ### Nan 0 to eliminate the 0/0.
  ### use naive situation
  #tmpnaive=matrix(rep((c-m$a)/(m$b-m$a),length(betahat)),nrow=k,ncol=n)
  tmp[is.nan(tmp)]= 0
  tmp
}

comp_postmean.tnormalmix = function(m,data){
  if(!is_normal(data$lik)){
    stop("Error: truncated normal mixture for non-normal likelihood is not yet implemented")
  }
  k=length(m$pi)
  n=length(data$x)
  A=1/(outer(m$sd^2,data$s^2,FUN="/")+1)
  B=1/outer(1/m$sd^2,1/data$s^2,FUN="+")
  ## try my_etruncnorm(1:3,2:4,matrix(0,3,4),matrix(1,3,4))
  result = my_etruncnorm(m$a,m$b,A*m$mean+t(t(1-A)*data$x),sqrt(B))
  ismissing = which(is.na(data$x) | is.na(data$s))
  if(length(ismissing)>0){result[,ismissing]=m$mean} 
  return(result)
}

comp_postsd.tnormalmix = function(m,data){
  if(!is_normal(data$lik)){
    stop("Error: truncated normal mixture for non-normal likelihood is not yet implemented")
  }
  k=length(m$pi)
  n=length(data$x)
  A=1/(outer(m$sd^2,data$s^2,FUN="/")+1)
  B=1/outer(1/m$sd^2,1/data$s^2,FUN="+")
  result = sqrt(my_vtruncnorm(m$a,m$b,t(A*m$mean+t(t(1-A)*data$x)),t(sqrt(B))))
  result = t(result)
  return(result)
}

comp_postmean2.tnormalmix = function(m,data){
  comp_postsd.tnormalmix(m,data)^2 + comp_postmean.tnormalmix(m,data)^2
}





