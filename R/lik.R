#' @title Likelihood object for normal error distribution
#' @description Creates a likelihood object for ash for use with normal error distribution
#' 
#' @examples 
#'    z = rnorm(100) + rnorm(100) # simulate some data with normal error
#'    ash(z,1,lik=normal_lik())
#' @export
normal_lik= function(){
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
#'    ash(z,1,lik=t_lik(df=4))
#' @export
t_lik = function(df){
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
#'    e = rnorm(1000) + log(rf(1000,df1=10,df2=10)) # simulate some data with log(F) error
#'    ash(e,1,lik=logF_lik(df1=10,df2=10))
#' @export
logF_lik = function(df1,df2){
  list(name = "logF",
       const = (length(unique(df1))==1) & (length(unique(df2))==1),
       lcdfFUN = function(x){plogf(x,df1=df1,df2=df2,log.p=TRUE)},
       lpdfFUN = function(x){dlogf(x,df1=df1,df2=df2,log=TRUE)})
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

