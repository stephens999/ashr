context("ashr with \"unit\" prior")

test_that("ash emits error if try to use prior=unit without VB", {
  expect_error(ash(rnorm(10),1,optmethod="mixEM",prior="unit"))
})

test_that("ash emits error if try to use prior<1 without VB", {
  expect_error(ash(rnorm(10),1,optmethod="mixEM",nullweight=0.5))
})
