test_that("compdens_conv for tnormal gives same results as compdens_conv for normal when a=-Inf, b=Inf", {
  gn = normalmix(c(0.5,0.5),c(0,0),c(0.1,1))
  gtn = tnormalmix(c(0.5,0.5),c(0,0),c(0.1,1),c(-Inf,-Inf),c(Inf,Inf))
  x=c(-10,2)
  s = c(1,2)
  data = list(x=x,s=s,lik=lik_normal(),alpha=0)
  expect_equal(comp_dens_conv(gn, data),comp_dens_conv(gtn,data))
})

test_that("log_compdens_conv for tnormal gives reasonable results", {
  gn = normalmix(c(0.5,0.5),c(0,0),c(0.1,1))
  gtn = tnormalmix(c(0.5,0.5),c(0,0),c(0.1,1),c(-Inf,-Inf),c(Inf,Inf))
  x=c(-10,2)
  s = c(1,2)
  data = list(x=x,s=s,lik=lik_normal(),alpha=0)
  expect_equal(exp(log_comp_dens_conv(gtn,data)),comp_dens_conv(gtn,data))
})
          
test_that("comp_cdf for tnormal gives same results as comp_cdf for normal when a=-Inf, b=Inf", {
  gn = normalmix(c(0.5,0.5),c(0,0),c(0.1,1))
  gtn = tnormalmix(c(0.5,0.5),c(0,0),c(0.1,1),c(-Inf,-Inf),c(Inf,Inf))
  y=c(-1,-5,0.5,1)
  expect_equal(comp_cdf(gn, y),comp_cdf(gtn,y))
})

test_that("compcdf_post for tnormal gives same results as compcdf_post for normal when a=-Inf, b=Inf", {
  gn = normalmix(c(0.5,0.5),c(0,0),c(0.1,1))
  gtn = tnormalmix(c(0.5,0.5),c(0,0),c(0.1,1),c(-Inf,-Inf),c(Inf,Inf))
  betahat = c(-1,-2,1,2)
  sebetahat = 1:4
  data = list(x=betahat,s=sebetahat,lik=lik_normal(),alpha=0)
  c = 0.5
  expect_equal(comp_cdf_post(gn,c,data),comp_cdf_post(gtn,c,data))
})

test_that("comp_postmean for tnormal gives same results as comp_postmean for normal when a=-Inf, b=Inf", {
  gn = normalmix(c(0.5,0.5),c(0,1),c(0.1,1))
  gtn = tnormalmix(c(0.5,0.5),c(0,1),c(0.1,1),c(-Inf,-Inf),c(Inf,Inf))
  betahat = c(-1,0.4,1.3,5)
  sebetahat = 1:4
  data = list(x=betahat,s=sebetahat,lik=lik_normal(),alpha=0)
  expect_equal(comp_postmean(gn,data),comp_postmean(gtn,data))
})

test_that("comp_postsd for tnormal gives same results as comp_postsd for normal when a=-Inf, b=Inf", {
  gn = normalmix(c(0.5,0.5),c(0,1),c(0.1,1))
  gtn = tnormalmix(c(0.5,0.5),c(0,1),c(0.1,1),c(-Inf,-Inf),c(Inf,Inf))
  betahat = c(-1,0.4,1.3,5)
  sebetahat = 1:4
  data = list(x=betahat,s=sebetahat,lik=lik_normal(),alpha=0)
  expect_equal(comp_postsd(gn,data),comp_postsd(gtn,data))
  expect_equal(comp_postmean2(gn,data),comp_postmean2(gtn,data))
})

# test_that("the my_etruncnorm in tnormal postmean calculation won't report error", {
#   beta = c(rep(0,100),rtruncnorm(100,a=5,b=10))
#   sebetahat = abs(rnorm(200,0,1))
#   betahat = rnorm(200,beta,sebetahat)
#   completeobs = (!is.na(betahat) & !is.na(sebetahat))
#   n=sum(completeobs)
#   
#   #Handling control variables
#   optmethod = 'mixIP'
#   randomstart=FALSE
#   gridmult=sqrt(2)
#   control=list()
#   control.default=list(K = 1, method=3, square=TRUE, step.min0=1, step.max0=1, mstep=4, kr=1, objfn.inc=1,tol=1.e-07, maxiter=5000, trace=FALSE)
#   controlinput = modifyList(control.default, control)
#   ##2. Generating mixture distribution
#   
#   mixsd = autoselect.mixsd(betahat[completeobs],sebetahat[completeobs],gridmult)
#   null.comp = which.min(mixsd) #which component is the "null"
#   
#   k = length(mixsd)
#   prior='uniform'
#   prior = setprior(prior,k,nullweight,null.comp)
#   pi = initpi(k,n,null.comp,randomstart)
#   
#   if(min(mixsd)>0){ #simply reflect the components
#     g = tnormalmix(c(pi,pi)/2,mean=0,sd=c(mixsd,mixsd),a=c(rep(-Inf,k),rep(0,k)),b=c(rep(0,k),rep(Inf,k)))
#     prior = rep(prior, 2)
#     pi = rep(pi, 2)
#   } else { #define two sets of components, but don't duplicate null component
#     null.comp=which.min(mixsd)
#     g = tnormalmix(c(pi,pi[-null.comp])/2,c(mixsd,mixsd[-null.comp]),c(rep(-Inf,k),rep(0,k-1)),c(rep(0,k),rep(Inf,k-1)))
#     prior = c(prior,prior[-null.comp])
#     pi = c(pi,pi[-null.comp])
#   }
#   m = g
#   A=outer(m$sd^2,sebetahat^2,FUN="/")/(outer(m$sd^2,sebetahat^2,FUN="/")+1)
#   B=1/outer(1/m$sd^2,1/sebetahat^2,FUN="+") 
#   alpha = (m$a-(t(t(A)*betahat)+(1-A)*m$mean))/sqrt(B)
#   beta =  (m$b-(t(t(A)*betahat)+(1-A)*m$mean))/sqrt(B)
#   #Flip the onese where both are positive, as the computations are more stable
#   #when both negative
#   flip = (alpha>0 & beta>0)
#   flip[is.na(flip)]=FALSE #deal with NAs
#   alpha[flip]= -alpha[flip]
#   beta[flip]=-beta[flip]
#   
#   #Fix a bug of quoting the truncnorm package
#   #E(X|a<X<b)=a when a==b as a natural result
#   #while etruncnorm would simply return NaN,causing PosteriorMean also NaN
#   isequal=(alpha==beta)
#   expect_equal(sum(is.na(isequal)),0)  
#   
# })


