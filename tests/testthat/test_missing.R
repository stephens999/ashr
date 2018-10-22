context("ashr with missing data")

test_that("missing data don't change results", {
  set.seed(11); z = rnorm(1000,0,2); s = rgamma(1000,10,10)
  z2 = c(z,rnorm(1000)); s2 = c(s, rep(Inf,1000))
  a = ash(z,s); a2 = ash(z2,s2)
  expect_equal(get_psd(a), get_psd(a2)[1:1000],tolerance=0.001)
  expect_equal(get_lfsr(a), get_lfsr(a2)[1:1000],tolerance=0.001)
})
