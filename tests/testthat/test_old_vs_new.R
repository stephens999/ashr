test_that("new results match old ones ", {
# these results were saved under v1.1.3 before introducing more stable
# calculation of likelihoods via first computing log-likelihood and normalizing
#  set.seed(1); z=rnorm(100); z.ash=ash(z,1); saveRDS(z.ash,file="tests/testthat/z.ash.test")
  set.seed(1); z=rnorm(100); z.ash=ash(z,1);
  oldres = readRDS("z.ash.test")
  expect_equal(modifyList(oldres,list(optmethod=NULL)),
               modifyList(z.ash,list(optmethod=NULL)))
})
