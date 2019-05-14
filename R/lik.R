#' @title Likelihood object for normal error distribution
#' @description Creates a likelihood object for ash for use with normal error distribution
#' 
#' @examples 
#'    z = rnorm(100) + rnorm(100) # simulate some data with normal error
#'    ash(z,1,lik=lik_normal())
#' @export
lik_normal = function(){
  list(name="normal",
       const = TRUE, #used to indicate whether the likelihood function is constant for all observations (some parts of ash only work in this case)
       lcdfFUN = function(x){stats::pnorm(x,log=TRUE)},
       lpdfFUN = function(x){stats::dnorm(x,log=TRUE)},
       etruncFUN = function(a,b){my_etruncnorm(a,b)},
       e2truncFUN = function(a,b){my_e2truncnorm(a,b)}
       )
}

#' @title Likelihood object for t error distribution
#' @description Creates a likelihood object for ash for use with t error distribution
#' @param df degree of freedom parameter of t distribution
#' 
#' @examples 
#'    z = rnorm(100) + rt(100,df=4) # simulate some data with t error
#'    ash(z,1,lik=lik_t(df=4))
#' @export
lik_t = function(df){
  list(name = "t",
      const = (length(unique(df))==1),
      lcdfFUN = function(x){stats::pt(x,df=df,log=TRUE)},
      lpdfFUN = function(x){stats::dt(x,df=df,log=TRUE)},
      etruncFUN = function(a,b){etrunct::e_trunct(a,b,df=df,r=1)},
      e2truncFUN = function(a,b){etrunct::e_trunct(a,b,df=df,r=2)}
      )
}

#' @title Likelihood object for logF error distribution
#' @description Creates a likelihood object for ash for use with logF error distribution
#' @param df1 first degree of freedom parameter of F distribution
#' @param df2 second degree of freedom parameter of F distribution
#' 
#' @examples 
#'    e = rnorm(100) + log(rf(100,df1=10,df2=10)) # simulate some data with log(F) error
#'    ash(e,1,lik=lik_logF(df1=10,df2=10))
#' @export
lik_logF = function(df1,df2){
  list(name = "logF",
       const = (length(unique(df1))==1) & (length(unique(df2))==1),
       lcdfFUN = function(x){plogf(x,df1=df1,df2=df2,log.p=TRUE)},
       lpdfFUN = function(x){dlogf(x,df1=df1,df2=df2,log=TRUE)},
       etruncFUN = function(a,b){
         matrix(mapply(my_etrunclogf,c(a),c(b),df1,df2),ncol=dim(a)[2])
         }
       )
}

#' @title Likelihood object for Poisson error distribution
#' @description Creates a likelihood object for ash for use with Poisson error distribution
#' @details Suppose we have Poisson observations \code{y} where \eqn{y_i\sim Poisson(c_i\lambda_i)}. 
#'    We either put an unimodal prior g on the (scaled) intensities \eqn{\lambda_i\sim g} 
#'    (by specifying \code{link="identity"}) or on the log intensities 
#'    \eqn{logit(\lambda_i)\sim g} (by specifying \code{link="log"}). Either way, 
#'    ASH with this Poisson likelihood function will compute the posterior mean of the 
#'    intensities \eqn{\lambda_i}.
#' @param y Poisson observations.
#' @param scale Scale factor for Poisson observations: y~Pois(scale*lambda).
#' @param link Link function. The "identity" link directly puts unimodal prior on Poisson
#'  intensities lambda, and "log" link puts unimodal prior on log(lambda).
#' 
#' @examples 
#'    beta = c(rnorm(100,50,5)) # prior mode: 50
#'    y = rpois(100,beta) # simulate Poisson observations
#'    ash(rep(0,length(y)),1,lik=lik_pois(y))
#'
#' @importFrom stats pgamma
#' @importFrom stats dgamma
#' 
#' @export
#'
lik_pois = function(y, scale=1, link=c("identity","log")){
  link = match.arg(link)
  if (link=="identity"){
    list(name = "pois",
         const = TRUE,
         lcdfFUN = function(x){pgamma(abs(x),shape=y+1,rate=scale,log.p=TRUE)-log(scale)},
         lpdfFUN = function(x){dgamma(abs(x),shape=y+1,rate=scale,log=TRUE)-log(scale)},
         etruncFUN = function(a,b){-my_etruncgamma(-b,-a,y+1,scale)},
         e2truncFUN = function(a,b){my_e2truncgamma(-b,-a,y+1,scale)},
         data=list(y=y, scale=scale, link=link))
  }else if (link=="log"){
    y1 = y+1e-5 # add pseudocount
    list(name = "pois",
         const = TRUE,
         lcdfFUN = function(x){pgamma(exp(-x),shape=y1,rate=scale,log.p=TRUE)-log(y1)},
         lpdfFUN = function(x){dgamma(exp(-x),shape=y1,rate=scale,log=TRUE)-log(y1)},
         etruncFUN = function(a,b){-my_etruncgamma(exp(-b),exp(-a),y1,scale)},
         e2truncFUN = function(a,b){my_e2truncgamma(exp(-b),exp(-a),y1,scale)},
         data=list(y=y, scale=scale, link=link))
  }
}

#' @title Likelihood object for Binomial error distribution
#' @description Creates a likelihood object for ash for use with Binomial error distribution
#' @details Suppose we have Binomial observations \code{y} where \eqn{y_i\sim Bin(n_i,p_i)}. 
#'    We either put an unimodal prior g on the success probabilities \eqn{p_i\sim g} (by specifying 
#'    \code{link="identity"}) or on the logit success probabilities \eqn{logit(p_i)\sim g} 
#'    (by specifying \code{link="logit"}). Either way, ASH with this Binomial likelihood function 
#'    will compute the posterior mean of the success probabilities \eqn{p_i}.
#' @param y Binomial observations
#' @param n Binomial number of trials
#' @param link Link function. The "identity" link directly puts unimodal prior on Binomial success
#'  probabilities p, and "logit" link puts unimodal prior on logit(p).
#' 
#' @examples 
#'    p = rbeta(100,2,2) # prior mode: 0.5
#'    n = rpois(100,10)
#'    y = rbinom(100,n,p) # simulate Binomial observations
#'    ash(rep(0,length(y)),1,lik=lik_binom(y,n))
#' @export
lik_binom = function(y,n,link=c("identity","logit")){
  link = match.arg(link)
  if (link=="identity"){
    list(name = "binom",
         const = TRUE,
         lcdfFUN = function(x){stats::pbeta(abs(x),shape1=y+1,shape2=n-y+1,log.p=TRUE)-log(n+1)},
         lpdfFUN = function(x){stats::dbeta(abs(x),shape1=y+1,shape2=n-y+1,log=TRUE)-log(n+1)},
         etruncFUN = function(a,b){-my_etruncbeta(-b,-a,y+1,n-y+1)},
         e2truncFUN = function(a,b){my_e2truncbeta(-b,-a,y+1,n-y+1)},
         data=list(y=y,n=n,link=link))
  }else if(link=="logit"){
    y1 = y+1e-5 # add pseudocount
    n1 = n+2e-5
    list(name = "binom",
         const = TRUE,
         lcdfFUN = function(x){stats::pbeta(1/(1+exp(-x)),shape1=y1,shape2=n1-y1,log.p=TRUE)+log(n1/(y1*(n1-y1)))},
         lpdfFUN = function(x){stats::dbeta(1/(1+exp(-x)),shape1=y1,shape2=n1-y1,log=TRUE)+log(n1/(y1*(n1-y1)))},
         etruncFUN = function(a,b){-my_etruncbeta(1/(1+exp(-b)),1/(1+exp(-a)),y1,n1-y1)},
         e2truncFUN = function(a,b){my_e2truncbeta(1/(1+exp(-b)),1/(1+exp(-a)),y1,n1-y1)},  
         data=list(y=y,n=n,link=link))
  }
}


#' @title Likelihood object for normal mixture error distribution
#' @description Creates a likelihood object for ash for use with normal mixture error distribution
#' @param pilik a k vector of mixture proportions (k is the number of mixture components), 
#'    or an n*k matrix that the j'th row the is mixture proportions for betahat_j
#' @param sdlik a k vector of component-wise standard deviations, 
#'    or an n*k matrix that the j'th row the is component-wise standard deviations for betahat_j
#' 
#' @examples 
#'    e = rnorm(100,0,0.8) 
#'    e[seq(1,100,by=2)] = rnorm(50,0,1.5) # generate e~0.5*N(0,0.8^2)+0.5*N(0,1.5^2)
#'    betahat = rnorm(100)+e
#'    ash(betahat, 1, lik=lik_normalmix(c(0.5,0.5),c(0.8,1.5)))
#' @export
lik_normalmix = function(pilik,sdlik){
  list(name="normalmix",
       const = (length(pilik)==length(sdlik)), #used to indicate whether the likelihood function is constant for all observations (some parts of ash only work in this case)
       lcdfFUN = function(x){
         cdfF = function(sdlik,x){stats::pnorm(x,mean=0,sd=sdlik)}
         if (length(pilik)==length(sdlik)){sdlik=matrix(sdlik,nrow=1)}
         sdlik = split(sdlik, rep(1:ncol(sdlik), each = nrow(sdlik)))
         log(Reduce("+",mapply('*', lapply(sdlik,cdfF,x=x), pilik, SIMPLIFY=FALSE)))},
       lpdfFUN = function(x){
         pdfF = function(sdlik,x){stats::dnorm(x,mean=0,sd=sdlik)}
         if (length(pilik)==length(sdlik)){sdlik=matrix(sdlik,nrow=1)}
         sdlik = split(sdlik, rep(1:ncol(sdlik), each = nrow(sdlik)))
         log(Reduce("+",mapply('*', lapply(sdlik,pdfF,x=x), pilik, SIMPLIFY=FALSE)))},
       etruncFUN = function(a,b){
         etrunc = function(sdlik,alpha,beta){my_etruncnorm(alpha,beta,mean=0,sd=sdlik)}
         if (length(pilik)==length(sdlik)){sdlik=matrix(sdlik,nrow=1)}
         sdlik = split(sdlik, rep(1:ncol(sdlik), each = nrow(sdlik)))
         Reduce("+",mapply('*', lapply(sdlik,etrunc,alpha=a,beta=b), pilik, SIMPLIFY=FALSE))},
       e2truncFUN = function(a,b){
         etrunc = function(sdlik,alpha,beta){my_e2truncnorm(alpha,beta,mean=0,sd=sdlik)}
         if (length(pilik)==length(sdlik)){sdlik=matrix(sdlik,nrow=1)}
         sdlik = split(sdlik, rep(1:ncol(sdlik), each = nrow(sdlik)))
         Reduce("+",mapply('*', lapply(sdlik,etrunc,alpha=a,beta=b), pilik, SIMPLIFY=FALSE))}
  )
}


# adds a default etruncFUN based on gen_etruncFUN, which uses numerical integration
add_etruncFUN = function(lik){
  if(is.null(lik$etruncFUN)){
    lik$etruncFUN = gen_etruncFUN(lik$lcdfFUN,lik$lpdfFUN)
  }
  if(is.null(lik$e2truncFUN)){ #for now add dummy function to return NA
    lik$e2truncFUN = function(a,b){tmp = rep(NA,length(a)); dim(tmp) = dim(a); return(tmp)} 
  }
  lik
}

is_normal = function(lik){lik$name=="normal"}
is_const = function(lik){lik$const}
get_name = function(lik){lik$name}

