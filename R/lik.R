normal_lik= function(){
  list(lcdfFUN = function(x){pnorm(x,log=TRUE)},
       lpdfFUN = function(x){dnorm(x,log=TRUE)},
       etruncFUN = function(a,b){my_etruncnorm(a,b)},
       e2truncFUN = function(a,b){my_e2truncnorm(a,b)}
       )
}

t_lik = function(df){
  list(lcdfFUN = function(x){ptgen(x,df=df,log=TRUE)},
      lpdfFUN = function(x){dtgen(x,df=df,log=TRUE)},
      etruncFUN = function(a,b){my_etrunct(a,b,df=df)},
      e2truncFUN = function(a,b){my_e2trunct(a,b,df=df)}
      )
}