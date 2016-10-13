test_that("some test about my_etrucnorm function", {
  expect_equal(-100,my_etruncnorm(-Inf,-100,0,1),tolerance=0.01)
  expect_equal(-100,my_etruncnorm(-Inf,-100,0,0))
  expect_equal(30,my_etruncnorm(30,100,0,0))
})
