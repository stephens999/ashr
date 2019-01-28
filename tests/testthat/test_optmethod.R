context("ashr optimization algorithms")

test_that("control is passed to optmethod correctly when method is mixIP", {
  testthat::skip_if_not_installed("REBayes")
  set.seed(1)
  z     <- rnorm(10,0,2) 
  z.ash <- ash(z,1,optmethod = "mixIP",control = list(rtol=1e-1),
                   outputlevel = 3)    
  expect_true(z.ash$fit_details$optreturn$control$rtol == 1e-1)
  
  z.ash <- ash(z,1,optmethod = "mixIP",outputlevel = 3)
  expect_true(z.ash$fit_details$optreturn$control$rtol == 1e-6)
})
