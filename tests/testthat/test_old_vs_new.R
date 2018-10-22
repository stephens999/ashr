context("ashr comparisons with previous versions")

test_that("new results match old ones ", {
# these results were saved under v1.1.3 before introducing more stable
# calculation of likelihoods via first computing log-likelihood and normalizing
# Then also under v1.1.10 after introducing svalue
# set.seed(1); z=rnorm(100,0,2); z.ash=ash(z,1); saveRDS(z.ash,file="tests/testthat/z.ash.test")
  set.seed(1); z=rnorm(100,0,2); z.ash=ash(z,1);
  oldres = readRDS("z.ash.test")
  expect_equal(get_pm(oldres),get_pm(z.ash),tolerance=0.001)
  expect_equal(get_psd(oldres),get_psd(z.ash),tolerance=0.001)
  expect_equal(get_lfsr(oldres),get_lfsr(z.ash),tolerance=0.001)
  expect_equal(get_fitted_g(oldres),get_fitted_g(z.ash),tolerance=0.001)
  expect_equal(get_logLR(oldres),get_logLR(z.ash),tolerance=0.001)
  expect_equal(get_loglik(oldres),get_loglik(z.ash),tolerance=0.001)
})
