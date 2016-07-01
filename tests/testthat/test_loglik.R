test_that("calc_null_loglik gives expected answer under normal", {
	set.seed(1); z=rnorm(100); s = rgamma(100,10,10); 
	alpha=0.5 
	expect_equal(sum(-alpha*log(s)+dnorm(z/(s^alpha),log=TRUE, sd=s^(1-alpha))),ashr::calc_null_loglik(list(x=z,s=s,v=NULL),alpha=alpha));
	expect_equal(sum(-log(s)+dnorm(z/s,log=TRUE, sd=1)),ashr::calc_null_loglik(list(x=z,s=s,v=NULL),alpha=alpha));
	
	alpha=0.1 
	expect_equal(sum(-alpha*log(s)+dnorm(z/(s^alpha),log=TRUE, sd=s^(1-alpha))),ashr::calc_null_loglik(list(x=z,s=s,v=NULL),alpha=alpha))
})


test_that("calc_null_loglik gives expected answer under t", {
        set.seed(1); z=rnorm(100); s = rgamma(100,10,10);
        alpha=0.5
        expect_equal(sum(-log(s) + dt(z/s,df=2,log=TRUE)),ashr::calc_null_loglik(list(x=z,s=s,v=2),alpha=alpha))
})

test_that("calc_logLR is 0 for when g is null", {
  set.seed(1); z=rnorm(100); 
  expect_equal(ashr::calc_logLR(g=unimix(1,0,0),list(x=z,s=1,v=NULL)),0)
  expect_equal(ashr::calc_logLR(g=unimix(1,0,0),list(x=z,s=1,v=2)),0)
})

test_that("calc_loglik returns warning when called with wrong model", {
  set.seed(1); z=rnorm(100); z.ash = ashr::ash(z,1)
  expect_warning(calc_loglik(z.ash,list(x=z,s=1,v=NULL),model="ET"))
})

test_that("logLR in ash object matches calc_logLR", {
  set.seed(1); z=rnorm(100); z.ash = ashr::ash(z,1)
  expect_equal(z.ash$logLR, calc_logLR(z.ash$fitted.g,list(x=z,s=1,v=NULL)))
  z.ash = ashr::ash(z,1,model="ET",df=3)
  expect_equal(z.ash$logLR, calc_logLR(z.ash$fitted.g,list(x=z,s=1,v=3),model="ET"))
})

test_that("sum of calc_vlogLR is same as calc_logLR", {
  set.seed(2); z=rnorm(100,0,2); z.ash = ashr::ash(z,1,df=4)
  expect_equal(sum(calc_vlogLR(z.ash$fitted.g,list(x=z,s=1,v=4))),z.ash$logLR)
})
