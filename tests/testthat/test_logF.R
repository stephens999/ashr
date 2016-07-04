test_that("logF error with df=10,10 gives similar answers to normal error", {
  set.seed(10)
  e = rnorm(100)+log(rf(100,df1=10,df2=10))
  e.ash = ash(e,1,lik=logF_lik(df1=10,df2=10))
  s = sd(log(rf(10000,10,10))) # get approximate standard deviation by simulation
  e.ash.norm = ash(e,s)
  expect_equal(e.ash.norm$PosteriorMean,e.ash$PosteriorMean,tol=0.05)
})