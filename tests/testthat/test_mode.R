context("ashr mode")

test_that("ash emits error if try to use mode and g", {
  expect_error(ash(rnorm(10),1,g=unimix(1,0,0),mode=0))
})

test_that("pm_on_zero gives expected results", {
  expect_equal(pm_on_zero(unimix(c(0.5,0.5),c(0,1),c(0,1))),c(TRUE,FALSE))
})


test_that("ash mode 1 gives same results as mode 0 but shifted by 1", {
  set.seed(1)
  z = rnorm(100,0,2)
  z.ash = ash(z,1)
  z.ash1 = ash(z+1,1,mode=1)
  expect_equal(get_pm(z.ash),get_pm(z.ash1)-1)
})
