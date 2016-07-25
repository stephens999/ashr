test_that("outputlevel works as expected", {
  set.seed(1); z=rnorm(10); z.ash=ash(z,1, outputlevel = c("fitted_g","logLR"))
  expect_null(z.ash$result)
  expect_type(z.ash$logLR,"double")
})
