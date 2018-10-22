context("ashr posterior sampling")

test_that("get_post_sample works as expected",{
  set.seed(1)
  n = 10 # number of observations
  se = 0.1
  nsamp = 1000 # number of samples per observation
  z = rnorm(n)
  z.ash = ash(z, se, "normal")
  samp = get_post_sample(z.ash, nsamp)
  # Check that the matrix of samples is of the correct dimensions:
  expect_equal(dim(samp), c(nsamp, n))
  # Check that the sampled posterior means are close to the true posterior means:
  expect_equal(colMeans(samp), z.ash$result$PosteriorMean, tolerance = 0.01)
  samp_sds = sqrt(apply(samp, 2, var))
  # Check that the sampled posterior SDs are close to the true posterior SDs:
  expect_equal(samp_sds, z.ash$result$PosteriorSD, tolerance = 0.01)
  
  u = runif(n)
  u.ash = ash(u, se, "uniform")
  u.samp = get_post_sample(u.ash, nsamp)
  expect_equal(dim(u.samp), c(nsamp, n))
  expect_equal(colMeans(u.samp), u.ash$result$PosteriorMean, tolerance = 0.01)
  u.samp_sds = sqrt(apply(u.samp, 2, var))
  expect_equal(u.samp_sds, u.ash$result$PosteriorSD, tolerance = 0.01)
})
