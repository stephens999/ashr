context("ashr with other likelihoods")

test_that("general likelihood with multiple df works", {
  df = c(rep(100,50),rep(2,50))
  s = rgamma(100,1,1)
  betahat = s*rt(n=100,df=df)
  data = set_data(betahat,s,lik = lik_normal(),alpha=0)
  expect_equal(is_normal(data$lik),TRUE)
  expect_equal(is_const(data$lik),TRUE)
  
  data =set_data(betahat,s,lik = lik_t(df),alpha=0)
  expect_equal(is_normal(data$lik),FALSE)
  expect_equal(is_const(data$lik),FALSE)
  
  data =set_data(betahat,s,lik = lik_t(1),alpha=0)
  expect_equal(is_normal(data$lik),FALSE)
  expect_equal(is_const(data$lik),TRUE)
})

test_that("general likelihood with multiple df works", {
  set.seed(10) 
  df = c(rep(100,50),rep(2,50))
  s = rgamma(100,1,1)
  betahat = s*rt(n=100,df=df)
  
  #calc_null_loglik(data)
  
  data =set_data(betahat,s,lik = lik_normal(),alpha=0)
  expect_equal(calc_null_loglik(data),sum(dnorm(betahat,sd=s,log=TRUE)))
  
  data =set_data(betahat,s,lik = lik_t(df),alpha=0)
  expect_equal(calc_null_loglik(data),sum(dt(betahat/s,df=df,log=TRUE)-log(s)))
  

})
