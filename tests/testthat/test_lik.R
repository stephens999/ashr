test_that("general likelihood with multiple df works", {
  set.seed(10) 
  df = c(rep(100,50),rep(2,50))
  s = rgamma(100,1,1)
  betahat = s*rt(n=100,df=df)
  
  #calc_null_loglik(data)
  
  data =set_data(betahat,s,lik = normal_lik(),alpha=0)
  expect_equal(calc_null_loglik(data),sum(dnorm(betahat,sd=s,log=TRUE)))
  
  data =set_data(betahat,s,lik = t_lik(df),alpha=0)
  expect_equal(calc_null_loglik(data),sum(dt(betahat/s,df=df,log=TRUE)-log(s)))
})