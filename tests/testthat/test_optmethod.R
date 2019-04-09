context("ashr optimization algorithms")

test_that("control is passed to optmethod correctly when method is mixIP", {

  # Since the Rmosek package on CRAN will not work with REBayes, here
  # I check whether the correct Rmosek package (the one downloaded
  # from mosek.com) is installed.
  testthat::skip_if_not_installed("REBayes")
  testthat::skip_if_not_installed("Rmosek")
  testthat::skip_if(is.element("mosek_attachbuilder",
                               getNamespaceExports("Rmosek")))
  set.seed(1)
  z     <- rnorm(10,0,2) 
  z.ash <- ash(z,1,optmethod = "mixIP",control = list(rtol=1e-1),
                   outputlevel = 3)    
  expect_true(z.ash$fit_details$optreturn$control$rtol == 1e-1)
  
  z.ash <- ash(z,1,optmethod = "mixIP",outputlevel = 3)
  expect_true(z.ash$fit_details$optreturn$control$rtol == 1e-6)
})
