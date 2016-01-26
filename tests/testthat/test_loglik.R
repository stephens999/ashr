test_that("calc_null_loglik gives expected answer under normal", {
	set.seed(1); z=rnorm(100); expect_equal(sum(log(dnorm(z))),ashr::calc_null_loglik(z,1,NULL))
})


test_that("calc_null_loglik gives expected answer under t", {
        set.seed(1); z=rnorm(100); expect_equal(sum(log(dt(z,df=2))),ashr::calc_null_loglik(z,1,df=2))
})

test_that("calc_logLR is 0 for when g is null", {
  set.seed(1); z=rnorm(100); expect_equal(ashr::calc_logLR(g=unimix(1,0,0),z,1,df=NULL),0)
  expect_equal(ashr::calc_logLR(g=unimix(1,0,0),z,1,df=2),0)
})

