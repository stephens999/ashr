test_that("optimization with weights matches expectations", {
  set.seed(1); z=rnorm(100,0,2) 
  z.ash= ash(z[1:50],1,optmethod="mixEM")
  z.ash.w = ash(z,1,optmethod="w_mixEM",weights = c(rep(1,50),rep(0,50)), g=get_fitted_g(z.ash))
  expect_equal(get_fitted_g(z.ash.w)$pi, get_fitted_g(z.ash)$pi, tol=0.001)
  
  testthat::skip_if_not_installed("REBayes")
  z.ash.w2 = ash(z,1,optmethod="mixIP",weights = c(rep(1,50),rep(0,50)), g=get_fitted_g(z.ash))
  expect_equal(get_fitted_g(z.ash.w2)$pi, get_fitted_g(z.ash)$pi, tol=0.001)
})
