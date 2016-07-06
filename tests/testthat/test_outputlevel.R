test_that("default output of posterior mean and sd occurs only when df=NULL", {
  set.seed(1); z=rnorm(10); 
  z.ash.n = ash(z,1)
  z.ash.t = ash(z,1,df=4)
  z.ash.t2 = ash(z,1,df=4,outputlevel=3)
  expect_true(!is.null(z.ash.n$PosteriorMean))
  expect_true(!is.null(z.ash.t2$PosteriorMean))
  expect_null(z.ash.t$PosteriorMean)
})