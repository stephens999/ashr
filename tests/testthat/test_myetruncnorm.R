context("my_etruncnorm")

test_that("my_etruncnorm returns expected results", {
  a    = c(-Inf,-Inf, 30,NA, 1, 1, 1,1,100,100,-100,-Inf,   0,-1,   1,  -2, -2)
  b    = c(-100,-100,100, 1,NA, 1, 1,1,Inf,Inf, -30, Inf, Inf, 1,   2,   2,  3)
  m    = c(   0,   0,  0, 1, 1,NA, 1,1,  0,  0,   0,   5,   0, 0,   0,-5e9,  0)
  sd   = c(   1,   0,  0, 1, 1, 1,NA,1,  1,  0,   0,   2,   1, 1,   1,   2,1e4)
  real = c(-100,-100, 30,NA,NA,NA,NA,1,100,100, -30,   5,0.80, 0,1.38,  -2,0.5)
  N = length(real)
  for (idx in 1:N){
    expect_equal(real[idx],my_etruncnorm(a[idx],b[idx],m[idx],sd[idx]),tolerance=0.01)
  }
  expect_equal(real,my_etruncnorm(a,b,m,sd),tolerance=0.01)
  real = matrix(real,N,4)
  m = matrix(m,N,4)
  sd = matrix(sd,N,4)
  expect_equal(real,my_etruncnorm(a,b,m,sd),tolerance=0.01)
  a=c(0,0)
  b=c(1,2)
  m = rbind(c(0,2,4),c(0,0,0))
  sd = 0
  real = rbind(c(0,1,1),c(0,0,0))
  expect_equal(real,my_etruncnorm(a,b,m,sd))
  
  expect_equal(my_etruncnorm(0,99,-2,3),truncnorm::etruncnorm(0,99,-2,3))
  expect_equal(my_etruncnorm(0,9999,-2,3),my_etruncnorm(0,Inf,-2,3),tol=1e-3)
  expect_error(my_etruncnorm(0, 1:2, mean = 0, sd = 1))
  expect_error(my_etruncnorm(1, 0, mean = 0, sd = 1))
  
  #TODO add test cases from pull request
})

context("my_vtruncnorm")

test_that("my_vtruncnorm returns expected results", {
  a    = c(-Inf,-Inf, 30,NA, 1, 1, 1,1,100,100,-100,-Inf,   0,  -1,   1,  -2,  -2)
  b    = c(-100,-100,100, 1,NA, 1, 1,1,Inf,Inf, -30, Inf, Inf,   1,   2,   2,   3)
  m    = c(   0,   0,  0, 1, 1,NA, 1,1,  0,  0,   0,   5,   0,   0,   0,-5e9,   0)
  sd   = c(   1,   0,  0, 1, 1, 1,NA,1,  1,  0,   0,   2,   2,   1,   1,   2, 1e4)
  real = c(   0,   0,  0,NA,NA,NA,NA,0,  0,  0,   0,   4,1.45,0.29,0.07,   0,2.08)
  N = length(real)
  for (idx in 1:N){
    expect_equal(real[idx],my_vtruncnorm(a[idx],b[idx],m[idx],sd[idx]),tolerance=0.01)
  }
  expect_equal(real, my_vtruncnorm(a, b, m, sd), tolerance = 0.01)
  real = matrix(real, N, 4)
  m = matrix(m, N, 4)
  sd = matrix(sd, N, 4)
  expect_equal(real, my_vtruncnorm(a, b, m, sd), tolerance = 0.01)
  a = c(0, 0)
  b = c(1, 2)
  m = rbind(c(0, 2, 4), c(0, 0, 0))
  sd = 0
  real = rbind(c(0, 0, 0), c(0, 0, 0))
  expect_equal(real, my_vtruncnorm(a,b,m,sd))
  
  expect_equal(my_vtruncnorm(-2, 3), truncnorm::vtruncnorm(-2, 3))
  expect_equal(my_vtruncnorm(6, 7, sd = 9), truncnorm::vtruncnorm(6, 7, sd = 9))
  expect_equal(my_vtruncnorm(0, 9999, -2, 3),
               my_vtruncnorm(0, Inf, -2, 3), tol = 1e-3)
})
