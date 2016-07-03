
#' @title my_e2trunct
#' @description Compute second moment of the truncated t. Uses results from O'Hagan, Biometrika, 1973
#' @param a left limit of distribution
#' @param b right limit of distribution
#' @param df degree of freedom of error distribution
#' @export
my_e2trunct = function(a,b,df,mean=0,sd=1){
  e_trunct(a,b,df,r=2)
}

#' @title my_etrunct
#' @description Compute second moment of the truncated t. Uses results from O'Hagan, Biometrika, 1973
#'
#' @param a left limit of distribution
#' @param b right limit of distribution
#' @param df degree of freedom of error distribution
#' @export
my_etrunct = function(a,b,df,mean=0,sd=1){
  e_trunct(a,b,df,r=1)
}


# These functions for computing moments of truncated t
# are based on results in O'Hagan, Biometrika 1973

ibeta = function(x,a,b){ pbeta(x,a,b)*beta(a,b) } #incomplete beta function

GG = function(r, v, q){
  p=v/(v+q^2)
  res = ifelse(q>0 | (r %% 2 != 0),
               ibeta(p,0.5*(v-r),0.5*(r+1)),
  2*beta(0.5*(v-r),0.5*(r+1)) - ibeta(p,0.5*(v-r),0.5*(r+1)))
  return(res)
}

# expectation of T^r when T has truncated t distribution on v df truncated at left by q
t2_1d = function(q,v,r){
  v^(0.5*r) * GG(r,v,q)/GG(0,v,q) 
}

# this function computes rth moment of t_v truncated at (a,b)
# works for v<r
e_trunct = function(a,b,v,r){
  mult=ifelse(a>0 & b>0,-1,1) # calculations more stable for negative
  #a and b, but in this case need to multiply result by (-1)^r
  aa = ifelse(a>0 & b>0, -b, a)
  bb = ifelse(a>0 & b>0, -a, b)
  
  mult^r * ifelse(aa==bb,aa^r,
         (pt(aa,v,lower.tail=FALSE)*t2_1d(aa,v,r) - 
              pt(bb,v,lower.tail=FALSE)*t2_1d(bb,v,r))/(pt(bb,v)-pt(aa,v)))
}



# e_trunct(3,4,4,1)
# my_etrunct(3,4,4)
# 
# e_trunct(-3,-4,4,1)
# my_etrunct(-3,-4,4)
# 
# e_trunct(3,4,4,2)
# my_e2trunct(3,4,4)
# 
# e_trunct(-3,-4,4,2)
# my_e2trunct(-3,-4,4)
# 
# e_trunct(-3,3,4,1)
# my_etrunct(-3,3,4)
# 
# e_trunct(-3,6,4,2)
# my_e2trunct(-3,6,4)
# 
# a = c(-3,-4,0)
# b = c(3,4,5)
# my_e2trunct(a,b,4)
# e_trunct(a,b,4,2)
