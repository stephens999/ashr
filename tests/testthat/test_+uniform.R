context("ashr with half-uniform mixture priors")

test_that("mixcompdist=+uniform gives all non-negative values for b and zero for a", {
  set.seed(1); z=rnorm(10); z.ash=ash(z,1,mixcompdist="+uniform")
  k = length(z.ash$fitted_g$pi)
  expect_true(all(z.ash$fitted_g$b >= rep(0,k)))
  expect_equal(z.ash$fitted_g$a,rep(0,k))
})

test_that("mixcompdist=-uniform gives all non-positive values for a and zero for b", {
  set.seed(1); z=rnorm(10); z.ash=ash(z,1,mixcompdist="-uniform")
  k = length(z.ash$fitted_g$pi)
  expect_equal(z.ash$fitted_g$b,rep(0,k))
  expect_true(all(z.ash$fitted_g$a <= 0))
})
