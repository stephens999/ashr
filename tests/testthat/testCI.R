context("ashr credible intervals")

test_that("CI result works for uniform prior", {
  set.seed(17)
  s = rgamma(100,10,10)
  z = rnorm(100,0,s+1)
  g = unimix(c(1),-1000,1000) #make prior very flat
  a = ash(z,s,g=g,fixg=TRUE)
  a2 = ash(z,s,df=4,g=g,fixg=TRUE)
  a3 = ash(z,s,alpha=1,g=g,fixg=TRUE)
 
  ci1=ashci(a,betaindex = 1:100)
  expect_equal(ci1[,1],z-1.96*s,tol=0.01)
  expect_equal(ci1[,2],z+1.96*s,tol=0.01)
  
  ci2 = ashci(a2,betaindex = 1:100)
  expect_equal(ci2[,1],z-qt(0.975,df=4)*s,tol=0.01)

  ci3 = ashci(a3,betaindex = 1:100)
  expect_equal(ci3[,1],z-1.96*s,tol=0.01)
})
