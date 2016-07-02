ibeta = function(x,a,b){ pbeta(x,a,b)*beta(a,b) } #incomplete beta function

GG = function(r, v, q){
  p=v/(v+q^2)
  if(q>0 | (r %% 2 != 0)){
    return(ibeta(p,0.5*(v-r),0.5*(r+1)))
  } else {
    return(2*beta(0.5*(v-r),0.5*(r+1)) - ibeta(p,0.5*(v-r),0.5*(r+1)))
  }
}

# expectation of T^r when T has truncated t distribution on v df truncated at left by q
t2_1d = function(q,v,r){
  v^(0.5*r) * GG(r,v,q)/GG(0,v,q) 
}

tt2 = function(a,b,v,r){
  ((1-pt(a,v))*t2_1d(a,v,r) - (1-pt(b,v))*t2_1d(b,v,r))/(pt(b,v)-pt(a,v))
}

tt2(3,4,4,2)
my_e2trunct(3,4,4)

tt2(3,4,4,1)
my_etrunct(3,4,4)

tt2(-3,3,4,1)
my_etrunct(-3,3,4)

tt2(-3,6,4,2)
my_e2trunct(-3,6,4)



