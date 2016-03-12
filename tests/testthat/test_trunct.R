test_that("my_e2trunct give sensible answers", {
  x = rt(1000000,df=4)
  expect_equal(mean(x[abs(x)<1]^2),my_e2trunct(-1,1,4),tolerance=0.01)
})

test_that("my_etrunct give sensible answers", {
  x = rt(1000000,df=4)
  expect_equal(mean(x[x>0 & x<2]),my_etrunct(0,2,4),tolerance=0.01)
  expect_equal(mean(x[x>1 & x<5]),my_etrunct(1,5,4),tolerance=0.01)
})

test_that("comp_postmean2.unimix gives right answer", {
  g= unimix(c(0.5,0.5),c(-1,0),c(1,2))
  comp_postmean2(g,c(1,2),c(1,2),4)
})