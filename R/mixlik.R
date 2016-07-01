
#' @title loglik_conv_mixlik
#'
#' @description find log likelihood of data, when the mixture m is convolved with l-comp normal mixture in betahat with mixture sd betahatsd, mixture proportion pilik.
#'
#' @param m mixture distribution
#' @param betahat an n vector
#' @param betahatsd an n-by-l matrix
#' @param v degree of freedom of error distribution
#' @param pilik an l vector
#' @param FUN default is "+"
#' @export
#'
loglik_conv_mixlik = function(m,betahat,betahatsd,v,pilik,FUN="+"){
  UseMethod("loglik_conv_mixlik")
}
#' @title loglik_conv.default
#'
#' @description The default version of \code{\link{loglik_conv}}.
#'
#' @param m mixture distribution
#' @param betahat an n vector
#' @param betahatsd an n-by-l matrix
#' @param v degree of freedom of error distribution
#' @param pilik an l vector
#' @param FUN default is "+"
#' @export
#'
loglik_conv_mixlik.default = function(m,betahat,betahatsd,v,pilik,FUN="+"){
  sum(log(dens_conv_mixlik(m,betahat,betahatsd,v,pilik,FUN)))
}

#compute the density of the components of the mixture m
#when convoluted with l-components normal mixture with standard deviation s
#or C (scale vector) multiplies scaled (se) student.t l-components mixture with df v
#with mixture proportion pilik
#the density is evaluated at x
#x is an n-vector
#s and pilik are n by l matrices
#v and c are l-vectors
#m is a mixture with k components
#output is a (k*l) by n matrix of densities

#' @title compdens_conv_mixlik
#' @description compute the density of the components of the mixture m
#'     when convoluted with l-components normal mixture with standard
#'     deviation s or C (scale vector) multiplies scaled (se)
#'     student.t l-components mixture with df v with mixture
#'     proportion pilik the density is evaluated at x
#'
#' @param m a mixture with k components
#' @param x an n vector
#' @param s normal mixture of sd(s), n-by-l matrix
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion , n-by-l matrix
#' @param FUN default is "+"
#' @return a (k*l) by n matrix of densities
#' @export
#'
#todo: C is missing for function input
compdens_conv_mixlik = function(m,x,s,v,pilik,FUN="+"){
  UseMethod("compdens_conv_mixlik")
}
compdens_conv_mixlik.default = function(m,x,s,v,pilik,FUN="+"){
  dens=NULL
  for (i in 1:dim(pilik)[2]){
    dens=rbind(dens,pilik[,i]*compdens_conv(m,x,s[,i],v[i],FUN))
  }
  return(dens)
}


#' @title dens_conv_mixlik
#' @description compute density of mixture m convoluted with
#'     l-components normal mixture of sd (s) or student t mixture with
#'     df v with mixture proportion pilik at locations x
#'
#' @param m mixture distribution
#' @param x an n vector
#' @param s normal mixture of sd(s), n-by-l matrix
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion , n-by-l matrix
#' @param FUN default is "+"
#' @export
#'
dens_conv_mixlik = function(m,x,s,v,pilik,FUN="+"){
  UseMethod("dens_conv_mixlik")
}
dens_conv_mixlik.default = function(m,x,s,v,pilik,FUN="+"){
  l=dim(pilik)[2]
  colSums(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik,FUN))
}


#' @title comppostprob_mixlik
#' @description compute the posterior prob that each observation came
#'     from each component of the mixture m,output a k by n vector of
#'     probabilities computed by weighting the component densities by
#'     pi and then normalizing,when likelihood is an l-components
#'     normal mixture or student t mixture with mixture proportion
#'     pilik
#'
#' @param m mixture distribution
#' @param x an n vector
#' @param s normal mixture of sd(s), n-by-l matrix
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion , n-by-l matrix
#' @export
comppostprob_mixlik=function(m,x,s,v,pilik){
  UseMethod("comppostprob_mixlik")
}
#' @export
comppostprob_mixlik.default = function(m,x,s,v,pilik){
  l=dim(pilik)[2]
  k=length(m$pi)
  tmp= (t(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik))/dens_conv_mixlik(m,x,s,v,pilik))
  ismissing = (is.na(x) | apply(is.na(s),1,sum))
  tmp[ismissing,]=rep(m$pi,l)/l
  group=rep(1:k,l)
  return(rowsum(t(tmp),group))
}


#' @title comppostprob_mixlik2
#' @description compute the posterior prob that each observation came
#'     from each component of the mixture m and the likelihood
#'     mixture, output a (k*l) by n vector of probabilities computed
#'     by weighting the component densities by pi and then
#'     normalizing, when likelihood is an l-components normal mixture
#'     or student t mixture with mixture proportion pilik.
#'
#' @param m mixture distribution
#' @param x an n vector
#' @param s normal mixture of sd(s), n-by-l matrix
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion , n-by-l matrix
#' @export
comppostprob_mixlik2=function(m,x,s,v,pilik){
  UseMethod("comppostprob_mixlik2")
}
#' @export
comppostprob_mixlik2.default = function(m,x,s,v,pilik){
  l=dim(pilik)[2]
  k=length(m$pi)
  tmp= (t(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik))/dens_conv_mixlik(m,x,s,v,pilik))
  ismissing = (is.na(x) | apply(is.na(s),1,sum))
  tmp[ismissing,]=rep(m$pi,l)/l
  return(t(tmp))
}



#' @title compcdf_post_mixlik
#' @description evaluate cdf of posterior distribution of beta at c
#' @param m prior on beta, a mixture
#' @param c location of evaluation
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @return it returns a (k*l) by n matrix
#' @export
compcdf_post_mixlik=function(m,c,betahat,sebetahat,v,pilik){
  UseMethod("compcdf_post_mixlik")
}
#' @export
compcdf_post_mixlik.default=function(m,c,betahat,sebetahat,v,pilik){
  cdf=NULL
  for (i in 1:dim(pilik)[2]){
    cdf=rbind(cdf,pilik[,i]*compcdf_post(m,c,betahat,sebetahat[,i],v[i]))
  }
  cdf
}

cdf_post_mixlik = function(m,c,betahat,sebetahat,v,pilik){
  UseMethod("cdf_post_mixlik")
}
cdf_post_mixlik.default=function(m,c,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik)*
            compcdf_post_mixlik(m,c,betahat,sebetahat,v,pilik))
}


#' @title postmean_mixlik
#' @description output posterior mean for beta for prior mixture m
#' @param m mixture distribution
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @export
#'
postmean_mixlik = function(m, betahat,sebetahat,v,pilik){
  UseMethod("postmean_mixlik")
}
postmean_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik) * comp_postmean_mixlik(m,betahat,sebetahat,v,pilik))
}


#' @title postmean2_mixlik
#' @description output posterior mean-squared value for beta for prior mixture m,given observations betahat, sebetahat, df v, from l-components mixture likelihood with mixture proportion pilik
#' @param m mixture distribution
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @export
postmean2_mixlik = function(m, betahat,sebetahat,v,pilik){
  UseMethod("postmean2_mixlik")
}
#' @export
postmean2_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik) * comp_postmean2_mixlik(m,betahat,sebetahat,v,pilik))
}

#' @title postsd_mixlik
#' @description output posterior sd for beta for prior mixture m
#' @param m mixture distribution
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @export
#'
postsd_mixlik = function(m,betahat,sebetahat,v,pilik){
  UseMethod("postsd_mixlik")
}
postsd_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  sqrt(postmean2_mixlik(m,betahat,sebetahat,v,pilik)-postmean_mixlik(m,betahat,sebetahat,v,pilik)^2)
}


#' @title comp_postmean2_mixlik
#' @description output posterior mean-squared value for beta
#' @param m mixture distribution
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @export
comp_postmean2_mixlik = function(m,betahat,sebetahat,v,pilik){
  UseMethod("comp_postmean2_mixlik")
}
#' @export
comp_postmean2_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  comp_postsd_mixlik(m,betahat,sebetahat,v,pilik)^2 +
    comp_postmean_mixlik(m,betahat,sebetahat,v,pilik)^2
}


#' @title comp_postmean_mixlik
#' @description output posterior mean for beta for each component of
#'     prior mixture m,given observations betahat, sebetahat, df v
#'     from l-components mixture likelihood with mixture proportion
#'     pilik
#' @param m mixture distribution
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @export
#'
comp_postmean_mixlik=function(m,betahat,sebetahat,v,pilik){
  UseMethod("comp_postmean_mixlik")
}
comp_postmean_mixlik.default=function(m,betahat,sebetahat,v,pilik){
  mean=NULL
  for (i in 1:dim(pilik)[2]){
    mean=rbind(mean,comp_postmean(m,betahat,sebetahat[,i],v[i]))
  }
  return(mean)
}


#' @title comp_postsd_mixlik
#'
#' @description output posterior sd for beta for each component of
#'     prior mixture m,given observations betahat, sebetahat, df v
#'     from l-components mixture likelihood with mixture proportion
#'     pilik
#' @param m mixture distribution
#' @param betahat the data
#' @param sebetahat the observed standard errors
#' @param v degree of freedom of error distribution
#' @param pilik mixture proportion
#' @export
#'
comp_postsd_mixlik=function(m,betahat,sebetahat,v,pilik){
  UseMethod("comp_postsd_mixlik")
}
comp_postsd_mixlik.default=function(m,betahat,sebetahat,v,pilik){
  sd=NULL
  for (i in 1:dim(pilik)[2]){
    sd=rbind(sd,comp_postsd(m,betahat,sebetahat[,i],v[i]))
  }
  return(sd)
}

