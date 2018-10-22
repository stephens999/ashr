context("truncated-normal computations")

test_that("gen_etruncFUN_single gives similar answers to etruncnorm and etrunct", {
  expect_equal(
    gen_etruncFUN_single(function(x){pnorm(x,log=TRUE)},
                          function(x){dnorm(x,log=TRUE)})(0,1),
    my_etruncnorm(0,1))
  expect_equal(
    gen_etruncFUN_single(function(x){pt(x,df=4,log=TRUE)},
                          function(x){dt(x,df=4,log=TRUE)})(0,1),
    my_etrunct(0,1,df=4))
})

test_that("gen_etruncFUN gives similar answers to etruncnorm and etrunct", {
  a=cbind(c(1,2),c(3,4))
  b=cbind(c(5,6),c(7,8))
  expect_equal(gen_etruncFUN(function(x){pt(x,df=4,log=TRUE)},
              function(x){dt(x,df=4,log=TRUE)})(a,b),
              my_etrunct(a,b,df=4))
})

test_that("ash with automatic truncfun gives similar answers to default", {
  set.seed(10); z=rnorm(10,0,4);
  a1 = ash(z,1)
  
  testlik = list(name="norm",lcdfFUN = function(x){pnorm(x,log=TRUE)},
       lpdfFUN = function(x){dnorm(x,log=TRUE)})
  a2 = ash(z,1,lik=testlik)
  expect_equal(a1$PosteriorMean,a2$PosteriorMean)
})
