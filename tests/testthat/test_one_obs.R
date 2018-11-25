context("ashr with 1 data sample")

test_that("ash works with one observation", {
  expect_error(ash(1,1),NA) #tests for no error
}
)
