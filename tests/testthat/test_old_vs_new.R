test_that("new results match old ones ", {
# these results were saved under v1.1.3 before introducing more stable
# calculation of likelihoods via first computing log-likelihood and normalizing
# Then also under v1.1.10 after introducing svalue
# set.seed(1); z=rnorm(100,0,2); z.ash=ash(z,1); saveRDS(z.ash,file="tests/testthat/z.ash.test")
  set.seed(1); z=rnorm(100,0,2); z.ash=ash(z,1);
  oldres = readRDS("z.ash.test")
  expect_equal(oldres$PosteriorMean,z.ash$res$PosteriorMean)
  expect_equal(oldres$PosteriorSD,z.ash$res$PosteriorSD)
  expect_equal(oldres$lfsr,z.ash$res$lfsr)
  expect_equal(oldres$fitted.g,z.ash$fitted.g)
  expect_equal(oldres$logLR,z.ash$logLR)
  expect_equal(oldres$loglik,z.ash$loglik)
})
