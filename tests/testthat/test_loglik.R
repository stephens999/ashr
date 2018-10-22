context("ashr likelihood computations")

test_that("calc_null_loglik gives expected answer under normal", {
	set.seed(1); z=rnorm(100); s = rgamma(100,10,10); 
	data1 = set_data(z,s,alpha=0.5)
	data2 = set_data(z,s,alpha=1)
	data3 = set_data(z,s,alpha=0.1)
        
	# likelihooods for different alpha should be same under null
	expect_equal(ashr::calc_null_loglik(data1),ashr::calc_null_loglik(data2))
	expect_equal(ashr::calc_null_loglik(data1),ashr::calc_null_loglik(data3))
	
	alpha=0.5
	expect_equal(sum(-alpha*log(s)+dnorm(z/(s^alpha),log=TRUE, sd=s^(1-alpha))),ashr::calc_null_loglik(data1));
	expect_equal(sum(-log(s)+dnorm(z/s,log=TRUE, sd=1)),ashr::calc_null_loglik(data1));
	
	alpha=0.1 
	expect_equal(sum(-alpha*log(s)+dnorm(z/(s^alpha),log=TRUE, sd=s^(1-alpha))),ashr::calc_null_loglik(data3))
})


test_that("calc_null_loglik gives expected answer under t", {
  set.seed(1); z=rnorm(100); s = rgamma(100,10,10);
  alpha=0.5
  data1 = set_data(z,s,lik_t(2),alpha=0.5)
  data2 = set_data(z,s,lik_t(2),alpha=1)
  expect_equal(ashr::calc_null_loglik(data1),ashr::calc_null_loglik(data2))
  expect_equal(sum(-log(s) + dt(z/s,df=2,log=TRUE)),ashr::calc_null_loglik(data1))
})

test_that("calc_logLR is 0 for when g is null", {
  set.seed(1); z=rnorm(100); 
  data1 = set_data(z,1)
  data2 = set_data(z,1,lik_t(2))
  expect_equal(ashr::calc_logLR(g=unimix(1,0,0),data1),0)
  expect_equal(ashr::calc_logLR(g=unimix(1,0,0),data2),0)
})

test_that("calc_loglik returns warning when called with wrong model", {
  set.seed(1); z=rnorm(100); z.ash = ashr::ash(z,1)
  data = set_data(z,1,alpha=1)
  expect_warning(calc_loglik(z.ash,data))
})

test_that("logLR in ash object matches calc_logLR", {
  set.seed(1); z=rnorm(100); z.ash = ashr::ash(z,1)
  data1 = set_data(z,1)
  expect_equal(z.ash$logLR, calc_logLR(z.ash$fitted_g,data1))
  z.ash = ashr::ash(z,1,alpha=1,df=3)
  data = set_data(z,1,lik_t(3),alpha=1)
  expect_equal(z.ash$logLR, calc_logLR(z.ash$fitted_g,data))
})

test_that("sum of calc_vlogLR is same as calc_logLR", {
  set.seed(2); z=rnorm(100,0,2); z.ash = ashr::ash(z,1,df=4)
  data = set_data(z,1,lik_t(4))
  expect_equal(sum(calc_vlogLR(z.ash$fitted_g,data)),z.ash$logLR)
})
