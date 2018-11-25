context("check ashr outputs")

test_that("outputlevel works as expected", {
  set.seed(1); z=rnorm(10); z.ash=ash(z,1, outputlevel = c("fitted_g","logLR"))
  expect_null(z.ash$result)
  expect_type(z.ash$logLR,"double")
})

test_that("penloglik produced correctly when fixg=TRUE and outputlevel=5", {
  set.seed(1); z=rnorm(10); z.ash=ash(z,1)
  g= get_fitted_g(z.ash)
  z.ash = ash(z,1,outputlevel = 5)
  z.ash2 = ash(z,1,fixg=TRUE,g=g,outputlevel = 5)
  expect_equal(z.ash$flash_data$penloglik,z.ash2$flash_data$penloglik)
})
