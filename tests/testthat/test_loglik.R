test_that("calc_null_loglik gives expected answer under normal", {
	set.seed(1); z=rnorm(100); s = rgamma(100,10,10); 
	alpha=0.5 
	expect_equal(sum(-alpha*log(s)+dnorm(z/(s^alpha),log=TRUE, sd=s^(1-alpha))),ashr::calc_null_loglik(z,s,NULL,alpha=alpha));
	expect_equal(sum(-log(s)+dnorm(z/s,log=TRUE, sd=1)),ashr::calc_null_loglik(z,s,NULL,alpha=alpha));
	
	alpha=0.1 
	expect_equal(sum(-alpha*log(s)+dnorm(z/(s^alpha),log=TRUE, sd=s^(1-alpha))),ashr::calc_null_loglik(z,s,NULL,alpha=alpha))
})


test_that("calc_null_loglik gives expected answer under t", {
        set.seed(1); z=rnorm(100); s = rgamma(100,10,10);
        alpha=0.5
        expect_equal(sum(-log(s) + dt(z/s,df=2,log=TRUE)),ashr::calc_null_loglik(z,s,df=2,alpha=alpha))
})

test_that("calc_logLR is 0 for when g is null", {
  set.seed(1); z=rnorm(100); expect_equal(ashr::calc_logLR(g=unimix(1,0,0),z,1,df=NULL),0)
  expect_equal(ashr::calc_logLR(g=unimix(1,0,0),z,1,df=2),0)
})

test_that("calc_loglik returns warning when called with wrong model", {
  set.seed(1); z=rnorm(100); z.ash = ashr::ash(z,1)
  expect_warning(calc_loglik(z.ash,z,1,NULL,model="ET"))
})

test_that("logLR in ash object matches calc_logLR", {
  set.seed(1); z=rnorm(100); z.ash = ashr::ash(z,1)
  expect_equal(z.ash$logLR, calc_logLR(z.ash$fitted.g,z,1,df=NULL))
  z.ash = ashr::ash(z,1,model="ET",df=2)
  expect_equal(z.ash$logLR, calc_logLR(z.ash$fitted.g,z,1,df=2,model="ET"))
})